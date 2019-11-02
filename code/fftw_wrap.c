#include "header.h"


static  fftw_plan       *plan, *iplan;
static  fftw_plan       *plan1D, *iplan1D;
static  fftw_complex    *work, *workT1, *workT2;
static  fftw_complex    *work1D;
static  fftw_plan        planT, iplanT;

static  int 		 N0;
static  int              M;
static  int 		 myid, np;
static  int              tmrFFT;


/* ---------------------------------------------------------------- */

void fft_init(geom_ptr geom,  fftw_complex *work_external){

/*
 * 	(nx*np)*ny  FFT computed in parallel (MPI FFTW)
 *	
 *	slab decomposed (along y) among numtasks MPI tasks
 *		
 *	ly  		         local ny
 *	lx			 local nx after transpose
 *	ly0			 local y start
 *	lx0			 local x start after transpose
 *	lsize			 total local size
 *		
 */

  int  ly, ly0, lx, lx0, lsize;      /* LOCAL DIMS */
  int  N, Ngrids, grid;


  /*-- get geometry information --*/ 

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  Ngrids = geom->Ngrids;
  N0     = geom->N0;
  N      = N0 * pow(2,Ngrids);

  work  = work_external;
  work1D  = (fftw_complex *) malloc( N * sizeof(fftw_complex));

 
  /*-- create plans --*/

  //fftw_init_threads();
  fftw_mpi_init();

  plan  = (fftw_plan *) malloc( Ngrids * sizeof(fftw_plan));
  iplan = (fftw_plan *) malloc( Ngrids * sizeof(fftw_plan));

  plan1D  = (fftw_plan *) malloc( Ngrids * sizeof(fftw_plan));
  iplan1D = (fftw_plan *) malloc( Ngrids * sizeof(fftw_plan));

  for (grid = 0; grid < Ngrids; grid++) {

    N = N0 * pow(2,grid);
 
    plan[grid] = fftw_mpi_plan_dft_2d(N, N, work, work, 
                     MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);

    iplan[grid] = fftw_mpi_plan_dft_2d(N, N, work, work, 
                     MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);

    plan1D[grid]  = fftw_plan_dft_1d(N, work1D, work1D, FFTW_FORWARD, FFTW_ESTIMATE);
    iplan1D[grid] = fftw_plan_dft_1d(N, work1D, work1D, FFTW_BACKWARD, FFTW_ESTIMATE);

  }

  /*-- print decomposition info --*/

  /*
   debug:
     fftwnd_mpi_local_sizes(plan, &ly, &ly0, &lx, &lx0, &lsize);
     printf("ID=%d   x:  %d->%d (%d),   y:  %d->%d (%d),  lsize= %d\n",
	  myid,  lx0, lx+lx0, lx,  ly0, ly+ly0, ly, lsize);

     fftwnd_mpi_destroy_plan(plan);
     fftwnd_mpi_destroy_plan(iplan); 
  */

  fft_extra_init(geom);

  tmrFFT = timer_set("FFT");

}

/* ---------------------------------------------------------------- */

void fft_wrap(fftw_complex *data, int grid, int dir){

  double    overN;
  long int  i, ntot;
  int       N;

  timer_on(tmrFFT);

  N = N0 * pow(2,grid);

  overN = 1./N;
  ntot = N*N/np;

  for (i=0; i<ntot; i++) {
    work[i][0] = data[i][0];
    work[i][1] = data[i][1];
  }


  if (dir == FORWARD)   fftw_execute(plan[grid]);

  if (dir == BACKWARD)  fftw_execute(iplan[grid]); 
 

  for (i=0; i<ntot; i++) {
    data[i][0] = work[i][0] * overN;
    data[i][1] = work[i][1] * overN;
  }


  timer_off(tmrFFT);

}


/* ---------------------------------------------------------------- */

void fft_wrap_1D(fftw_complex *data, int grid, int dir){

  double    overN;
  long int  i;
  int       N;

  timer_on(tmrFFT);

  N = N0 * pow(2,grid);

  overN = sqrt(1./N);

  for (i=0; i<N; i++) {
    work1D[i][0] = data[i][0];
    work1D[i][1] = data[i][1];
  }


  if (dir == FORWARD)   fftw_execute(plan1D[grid]);

  if (dir == BACKWARD)  fftw_execute(iplan1D[grid]); 
 

  for (i=0; i<N; i++) {
    data[i][0] = work1D[i][0] * overN;
    data[i][1] = work1D[i][1] * overN;
  }


  timer_off(tmrFFT);

}

/* ---------------------------------------------------------------- */

void fft_time_init(post_ptr post){

  /*-- get geometry information --*/ 

  M = post->time_slices;

  /*-- allocate work array --*/ 

  workT1  = (fftw_complex *) malloc( M * sizeof(fftw_complex));
  workT2  = (fftw_complex *) malloc( M * sizeof(fftw_complex));

  /*-- create plans --*/

  //planT  = fftw_create_plan(M, FFTW_FORWARD, FFTW_ESTIMATE);
  //iplanT = fftw_create_plan(M, FFTW_BACKWARD, FFTW_ESTIMATE);


}


/* ---------------------------------------------------------------- */

void fft_time_wrap(fftw_complex *data, int dir){

  double    overM;
  long int  i, ntot, m;

  timer_on(tmrFFT);

  overM = 1./sqrt(M);
  ntot = N0*N0/np;

  for (i=0; i<ntot; i++) {

     for (m=0; m<M; m++)  {
       workT1[m][0] = data[m*ntot + i][0];
       workT1[m][1] = data[m*ntot + i][1];
     }

     //if (dir == FORWARD)  fftw_one(planT, workT1, workT2);
     //if (dir == BACKWARD) fftw_one(iplanT, workT1, workT2);

     for (m=0; m<M; m++)  {
        data[m*ntot + i][0] = workT2[m][0] *overM;
        data[m*ntot + i][1] = workT2[m][1] *overM;
     }

  }

  timer_off(tmrFFT);

}

/* ---------------------------------------------------------------- */
