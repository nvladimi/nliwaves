#include "header.h"

static  fftw_complex   *psi;           // main array for main variable
static  fftw_complex   *psihat;        // fft of psi
static  fftw_complex   *work;          // work array for fft and regridding
static  double         *dt;            // timesteps for all grids
static  int            *n;             // number of steps at each grid
static  char            msg[80];       // misc messages written to file
static  int             Ngrids;        // number of grids



static  void   main_allocate_arrays(int N0);

static  void   main_init_grids  (geom_ptr geom, ctrl_ptr ctrl);

static  int    is_coarsest_step();
static  double addtime_one_step(int grid);


/*------------------------------------------------------------------------*/


int main(int argc, char **argv)
{
  ctrl_str      ctrl;
  geom_str      geom;
  phys_str      phys;
  diag_str      diag;
  ic_str        ic;

  int           dnInfo, dnData, dnDiag;        // frequency of output
  int           njig;                          // frequency on vacuum layers
  int           noff;                          // width of vacuum layers
  int           tmrInit, tmrComp, tmrPost;     // timers
  int           tmrFLAG;                       // flag to output timers

  int           grid;                          // current grid number
  int           nio;                           // current IO number
  double        t;                             // current time
  double        tmax;                          // finish time

  int           collapses;
  int           diag_averaging;

  int           required, provided;

  /*-- read parameters and store them in structures --*/

  MPI_Init(&argc, &argv);
  //MPI_Init_thread(&argc, &argv, required, &provided);


  parameters(&argc, &argv, &geom, &ctrl, &phys, &diag);


  /*-- set grids and timesteps --*/

  main_init_grids(&geom, &ctrl);
  grid =  geom.grid; 


  /*-- set current time and IO frequencies --*/

  t       =  n[0]*dt[0];
  tmax    =  ctrl.tmax;
  dnInfo  =  ctrl.dnInfo;
  dnDiag  =  ctrl.dnDiag;
  dnData  =  ctrl.dnData;
  njig    =  ctrl.njig;
  noff    =  ctrl.noff;

  if (dnDiag < 0) {
    dnDiag = -dnDiag;
    diag_averaging = 1;
  } else diag_averaging = 0;

  /*-- set timers --*/

  timers_init(ctrl.runname, t);

  tmrInit = timer_set("initialization");
  tmrComp = timer_set("computing");
  tmrPost = timer_set("postproc");


  /*-- initialization --*/

  timer_on(tmrInit);

  main_allocate_arrays(geom.N0);

  io_init        (&ctrl, &geom);
  grids_init     (&geom, work);
  fft_init       (&geom, work);
  evolve_init    (&geom, &phys, psi, psihat);

  info_init      (&ctrl, &geom, &diag, psi, psihat);
  diag_init      (&ctrl, &geom, &diag, &phys, psi, psihat);


#ifdef FIBER_WALL
  mask_init      (&geom, &phys, psi);
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  sprintf(msg, "#msg: all modules are initialized\n");
  io_message(msg,0);
 
  collapses = collapses_init(&ctrl, &geom, &diag, &phys);

  nio = io_read_data(psi, geom.grid);
  sprintf(msg, "#msg: data %04d restored at t = %18.12e\n", nio, t);
  io_message(msg,0);

  diag_ic(psi, psihat, grid, nio, t);

  timer_off(tmrInit);
  timers_out(t);
  tmrFLAG = 0;

  sprintf(msg, "#msg: starting simulation\n\n");
  io_message(msg,0);


  /*-- main loop --*/

  while (t<tmax) {

    /*-- if jigging, make one linear step --*/

    //sprintf(msg, "#test: n=%d  ", n[0]);
    //io_message(msg,0);

    if ( (njig>0) && (is_coarsest_step()) && (n[0]>0) && (n[0]%njig==0) ) {

      //sprintf(msg, "#test: linear step\n");
      //io_message(msg,0);

      evolve_biglin(grid,  noff*dt[0]);
      n[0] += noff;
      t = dt[0]*n[0];
      continue;
    }

    //sprintf(msg, "#test: full step\n");
    //io_message(msg,0);



    /*-- computing --*/

    timer_on(tmrComp);

    grid = spectrum_regrid(grid, n[grid], t);

    if (grid == -1) {
      sprintf(msg, "#msg: additional finer grids needed, t = %18.12e\n", t);
      io_message(msg,0);
      info_output(grid, t);
      MPI_Finalize();
      exit(0); 
    }

#ifdef FIBER_WALL
    mask_apply(grid, dt[grid]);
#endif

    evolve_one_step(grid, dt[grid]);

    t = addtime_one_step(grid);

    timer_off(tmrComp);


    /*-- post-processing --*/

    timer_on(tmrPost);

    if (collapses) {
      collapses_prepare(psi, grid);
      collapses_update(t);
    }


    if ( is_coarsest_step() ) {

      if ( n[0]%dnInfo == 0 ) {

	info_output(grid, t);
	if (collapses){
	  //collapses_prepare(psi, grid);
	  collapses_addnew(t);
	  collapses_update(t);
	  collapses_remove(t);
	}

      }

      if ( n[0]%(dnInfo*dnDiag) == 0 )  {
	diag_compute(psi, psihat, grid);
	if (!diag_averaging) diag_output(nio+1);
      }

      if ( n[0]%(dnInfo*dnData) == 0 ) {
	nio = io_save_data(psi, grid);
        if (diag_averaging) diag_output(nio);
	io_save_tag(nio, n[0], grid, t);
	tmrFLAG = 1;
      }
    } 



    timer_off(tmrPost);


    /*-- output timers, if it was data step --*/ 

    if (tmrFLAG) {timers_out(t); tmrFLAG = 0;}

  }

  /*-- end of main loop --*/


  MPI_Finalize();
  exit(0); 
}

/*------------------------------------------------------------------------*/


void main_allocate_arrays(int N0){

  int N, np, myid;

  N = N0 * pow(2, Ngrids-1);

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  psi    = (fftw_complex *) malloc( N*N/np * sizeof(fftw_complex));
  psihat = (fftw_complex *) malloc( N*N/np * sizeof(fftw_complex));
  work   = (fftw_complex *) malloc( N*N/np * sizeof(fftw_complex));


}


/*------------------------------------------------------------------------*/

void main_init_grids(geom_ptr geom, ctrl_ptr ctrl){

  double dx0;
  int g;

  Ngrids = geom->Ngrids;  

  n  = (int*) malloc(Ngrids * sizeof(int));
  dt = (double*) malloc(Ngrids * sizeof(double));

  dx0     =  geom->Lx/geom->N0; 

  if (ctrl->dt0) dt[0] = ctrl->dt0;
  else           dt[0] = ctrl->CFL * dx0*dx0;

  n[0]    =  ctrl->n0;

  for (g=1; g<Ngrids; g++) {
    dt[g] = dt[g-1]/4; 
    n[g]=0;
  }

}


/*------------------------------------------------------------------------*/

int is_coarsest_step(){

  int g;

  for (g=1; g<Ngrids; g++)  if (n[g]) return(0);

  return(1);

}

/*------------------------------------------------------------------------*/

double addtime_one_step(int grid){

  int    g;
  double t;

  n[grid]++;

  t = 0;  

  for (g=grid; g>0; g--) {

    if (n[g] == 4) {
      n[g] = 0;
      n[g-1] ++;
    }

    t += dt[g]*n[g];

  }

  t += dt[0]*n[0];

  //sprintf(msg, "#dbg:  t = %18.12e,  n = %4d %d\n", t, n[0], n[1]);
  //io_message(msg,0);
 
  return(t);
  
}

/*------------------------------------------------------------------------*/
