#include "header.h"

int main(int argc, char **argv)
{
  fftw_complex 	**A, **D, **work;
  fftw_complex  *df0, *df1, *df2;

  ctrl_str  ctrl;
  geom_str  geom;
  ic_str    ic;

  int grid = 0;
  int gp   = 0;

  double   L = 1;
  int      N;
  int      myid, np;
  int      nio;

  char   filename[80];
  FILE   *thefile;
  int     i;

  /* --------------------------------------------- */

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if(argc != 3) {
    if (myid==0) printf("\n\t usage: %s run_name  N_grid \n\n", argv[0]);
    MPI_Finalize();
    exit(1);
  }

  /* set parameters of test field*/

  memset(&geom, 0, sizeof(&geom));
  memset(&ctrl, 0, sizeof(&ctrl));
  memset(&ic,   0, sizeof(&ic));

  strcpy(ctrl.runname, argv[1]);
  N = atoi(argv[2]);

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


  /* major array allocation */

  arr2D_init(&geom, gp);

  A       = arr2D_create(grid, 0);
  D       = arr2D_create(grid, 0);
  work    = arr2D_create(grid, 0);


  /* initialization */

  fft_init(&geom, work[0]);
  io_init(&ctrl, &geom);
  ic_set(A[0], &geom, &ic);

/* ------------------------------------------------- */

  /* field variable f(x,y) */

  //nio = io_save_data(A[0], grid);
  //if (myid==0)  printf("%4d:  f(x,y) \n", nio);


  /* fft transform forth and back */

  arr2D_copy(A, D, grid);

  for (i=0; i<10; i++)  fft_test(D[0], grid);

  //nio = io_save_data(D[0], grid);
  //if (myid==0)  printf("%4d:  f(x,y)    FFT fransformed twice - 0 cycles\n", nio);
  if (myid==0)  printf("FFT fransformed twice - 10 cycles, no binary output\n");

  MPI_Finalize();
  exit(0); 
}
