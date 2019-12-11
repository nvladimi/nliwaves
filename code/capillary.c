#include "header.h"

static  int 		 myid, np;
static  int 		 N0;

static  fftw_complex  *Psi;


/* --------------------------------------------- */


int main(int argc, char **argv)
{
  fftw_complex 	**work;
 
  ctrl_str  ctrl;
  geom_str  geom;
  phys_str  phys;
  ic_str    ic;

  int grid = 0;

  double   L = 1;
  int      N = 128;
  int      Nx, Ny;
  int      nio;

  char   filename[80];
  FILE   *thefile;
  int     i;

  /* --------------------------------------------- */

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if(argc != 2) {
    if (myid==0) printf("\n\t usage: %s run_name\n\n", argv[0]);
    MPI_Finalize();
    exit(1);
  }

  /* set parameters of test field*/

  memset(&geom, 0, sizeof(&geom));
  memset(&ctrl, 0, sizeof(&ctrl));
  memset(&ctrl, 0, sizeof(&phys));
  memset(&ic,   0, sizeof(&ic));

  strcpy(ctrl.runname, argv[1]);

  strcpy(ic.type, "ellipse");

  ctrl.restart = 0;

  ic.BumpHeight = 1;
  ic.BumpX = L/2;
  ic.BumpY = L/2;
  ic.BumpA = L/6;
  ic.BumpB = L/8;

  geom.Ngrids = 1;
  geom.grid = 0;
  geom.N0 = N;
  geom.Lx = L;
  geom.Ly = L;

  phys.coefA = 1;
  phys.coefB = 1;

  
  N0 = N;

  /* major array allocation */

  arr2D_init(&geom, 0);

  work    = arr2D_create(grid, 0);

  /*--- create suppementary arrays ---*/

  Nx = N0*pow(2, geom.Ngrids-1);
  Ny = Nx/np;

  Psi = (fftw_complex *) malloc( Nx*Ny * sizeof(fftw_complex));

  /* initialization */

  fft_init(&geom, work[0]);
  io_init(&ctrl, &geom);
  ic_set(Psi, &geom, &ic);

  rhs_init(&geom, &phys, Psi);
  rhs_init_S(Psi);


/* ------------------------------------------------- */

  /* field variable f(x,y) */

  nio = io_save_data(Psi, grid);
  if (myid==0)  printf("%4d:  capillary Psi\n", nio);

  rhs_compute(grid);

  nio = io_save_data(Psi, grid);
  if (myid==0)  printf("%4d:  capillary dPsi/dt\n", nio);

  MPI_Finalize();
  exit(0); 
}

/* ---------------------------------------------------------------- */

