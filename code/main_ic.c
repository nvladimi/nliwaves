#include "header.h"

static int    np, myid;
       char   msg[80];

void ic_header(geom_ptr geom, ctrl_ptr ctrl, ic_ptr ic);

void ic_parameters(int *argc, char ***argv, 
		geom_ptr geom, ctrl_ptr ctrl, ic_ptr ic);

/*------------------------------------------------------------------------*/


int main(int argc, char **argv)
{
  fftw_complex   *psi;           // main array for main variable
  fftw_complex   *work;          // main array for main variable

  ctrl_str      ctrl;
  geom_str      geom;
  phys_str      phys;
  diag_str      diag;
  ic_str        ic;

  double        t = 0;
  int           grid = 0;        // current grid number
  int           nio = 0;         // current IO number

  int           N0;

  /*-- read parameters and store them in structures --*/

  MPI_Init(&argc, &argv);

  ic_parameters(&argc, &argv, &geom, &ctrl, &ic);

  N0 = geom.N0;

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  psi    = (fftw_complex *) malloc( N0*N0/np * sizeof(fftw_complex));
  work   = (fftw_complex *) malloc( N0*N0/np * sizeof(fftw_complex));


  io_init(&ctrl, &geom);

  fft_init(&geom, work);

  ic_set(psi, &geom, &ic);

  nio = io_save_data(psi, grid);

  io_save_tag(nio, nio, grid, t);


  MPI_Finalize();
  exit(0); 

}

/*------------------------------------------------------------------------*/



void ic_parameters(int *argc, char ***argv, 
		geom_ptr geom, ctrl_ptr ctrl, ic_ptr ic)
{
  char       line[200], param[200], value[200];
  char       filename[80];  
  FILE      *thefile;
  double     pi;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if(argc[0] != 2) {
    if(!myid){printf("\n\t usage: %s run_name\n\n", argv[0][0]);}
    MPI_Finalize();
    exit(1);
  }

  /*-- zero out structures --*/

  memset(geom, 0, sizeof(*geom));
  memset(ctrl, 0, sizeof(*ctrl));
  memset(ic,   0, sizeof(*ic));

  /*-- get current run name and the name of parameter file --*/

  strcpy(ctrl->runname, argv[0][1]);
  strcpy(filename, ctrl->runname);
  strcat(filename, ".ic");

 
  /*-- master processor reads the input file and broadcasts the parameters --*/

  if (myid==0){

    /*-- reading input file --*/

    thefile = fopen (filename, "rt");

    if (thefile == NULL) {
      printf ("\n  Cannot open \"%s\".\n\n", filename);
      MPI_Finalize();
      exit(1);
    }
 
    while(fgets(line, 200, thefile) != NULL){

      sscanf(line, "%s  %s", param, value);

      if( strcmp(param,"N0") == 0 )	  geom->N0            = atoi(value);
      if( strcmp(param,"Lx") == 0 )	  geom->Lx            = atof(value);
      if( strcmp(param,"Ly") == 0 )	  geom->Ly            = atof(value);

      if( strcmp(param,"ic_type") == 0 )  strcpy(ic->type,value);

      if( strcmp(param,"ic_n")  == 0 )    ic->NumberOfBumps   = atoi(value);
      if( strcmp(param,"ic_h")  == 0 )	  ic->BumpHeight      = atof(value);
      if( strcmp(param,"ic_r")  == 0 )	  ic->BumpRadius      = atof(value);
      if( strcmp(param,"ic_x0") == 0 )	  ic->BumpX           = atof(value);
      if( strcmp(param,"ic_y0") == 0 )    ic->BumpY           = atof(value);
      if( strcmp(param,"ic_a")  == 0 )	  ic->BumpA           = atof(value);
      if( strcmp(param,"ic_b")  == 0 )    ic->BumpB           = atof(value);
      if( strcmp(param,"ic_chirp") == 0)  ic->Chirp           = atof(value);
      if( strcmp(param,"ic_seed")  == 0 ) ic->Seed            = atoi(value);
      if( strcmp(param,"ic_kmin") == 0 )  ic->kmin            = atof(value);
      if( strcmp(param,"ic_kmax") == 0 )  ic->kmax            = atof(value);
      if( strcmp(param,"ic_npart") == 0 ) ic->numParticles    = atof(value);
      if( strcmp(param,"ic_ncond") == 0 ) ic->nCondensate     = atof(value);


    }

    fclose(thefile);

  }

    /*-- derived parameters --*/

    geom->Ngrids = 1;
    geom->grid = 0;

    pi =  acos(0)*2;
 
    if (geom->Lx < 0) geom->Lx = -pi*geom->Lx;
    if (geom->Ly < 0) geom->Ly = -pi*geom->Ly;

    ic_header(geom, ctrl, ic);


  /*-- broadcast the parametrers --*/

  MPI_Bcast( geom, sizeof(geom_str), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast( geom, sizeof(ctrl_str), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(   ic, sizeof(ic_str),   MPI_BYTE, 0, MPI_COMM_WORLD);

}



/*-----------------------------------------------------------*/

void ic_header(geom_ptr geom, ctrl_ptr ctrl, ic_ptr ic)
{

  char       filename[80];  
  FILE      *thefile;
 
  if (myid != 0) return;


  strcpy (filename, ctrl->runname);
  strcat (filename, ".out");  thefile = fopen (filename, "a");

  fprintf(thefile, "\n");
  fprintf(thefile, "# SETTING INITIAL CONDITIONS\n\n");  
  fprintf(thefile, "# Processes: %d\n", np);  
  fprintf(thefile, "# Run name:  %s\n", ctrl->runname );

  fprintf(thefile, "#\n# Domain geometry: \n");
  fprintf(thefile, "#     Ngrids= %d\n",  geom->Ngrids);
  fprintf(thefile, "#     grid  = %d\n",  geom->grid);
  fprintf(thefile, "#     N0    = %d\n",  geom->N0);
  fprintf(thefile, "#     Lx    = %f\n",  geom->Lx);
  fprintf(thefile, "#     Ly    = %f\n",  geom->Ly);


  fprintf(thefile, "#\n# Initial conditions: %s\n",    ic->type);
  fprintf(thefile, "#     number of bumps  n  = %d\n", ic->NumberOfBumps);
  fprintf(thefile, "#     bump height      h  = %f\n", ic->BumpHeight);
  fprintf(thefile, "#     bump radius      r  = %f\n", ic->BumpRadius);
  fprintf(thefile, "#     bump x-coord     x0 = %f\n", ic->BumpX);
  fprintf(thefile, "#     bump y-coord     y0 = %f\n", ic->BumpY);
  fprintf(thefile, "#     bump radius A    a  = %f\n", ic->BumpA);
  fprintf(thefile, "#     bump radius B    b  = %f\n", ic->BumpB);
  fprintf(thefile, "#     random seed    seed = %d\n", ic->Seed);
  fprintf(thefile, "#     spectrum       kmin = %f\n", ic->kmin);
  fprintf(thefile, "#     spectrum       kmax = %f\n", ic->kmax);
  fprintf(thefile, "#     spectrum      npart = %f\n", ic->numParticles);
  fprintf(thefile, "#     spectrum      ncond = %f\n", ic->nCondensate);


  fprintf(thefile, "\n");

  fclose(thefile);

}

/*-----------------------------------------------------------*/








