#include "header.h"

int main(int argc, char **argv)
{
  fftw_complex 	**A, **Agp, **D, **work;
  fftw_complex  *df0, *df1, *df2;

  ctrl_str  ctrl;
  geom_str  geom;
  ic_str    ic;

  int grid = 0;
  int gp   = 2;

  double   L = 1;
  int      N = 128;
  int      myid, np;
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


  /* major array allocation */

  arr2D_init(&geom, gp);

  A       = arr2D_create(grid, 0);
  Agp     = arr2D_create(grid, 2);
  D       = arr2D_create(grid, 0);
  work    = arr2D_create(grid, 0);

  df0  = (fftw_complex *) malloc( N * sizeof(fftw_complex));
  df1  = (fftw_complex *) malloc( N * sizeof(fftw_complex));
  df2  = (fftw_complex *) malloc( N * sizeof(fftw_complex));


  /* initialization */

  fft_init(&geom, work[0]);
  io_init(&ctrl, &geom);
  ic_set(A[0], &geom, &ic);

/* ------------------------------------------------- */

  /* field variable f(x,y) */

  nio = io_save_data(A[0], grid);
  if (myid==0)  printf("%4d:  f(x,y) \n", nio);


  /* x-derivative */

  arr2D_copy_gp(A, Agp, grid);
  arr2D_deriv_x(Agp, D, grid);
  nio = io_save_data(D[0], grid);
  if (myid==0)  printf("%4d:  df/dx     (central differences)\n", nio);

  /* y-derivative */

  arr2D_copy_gp(A, Agp, grid);
  arr2D_deriv_y(Agp, D, grid);
  nio = io_save_data(D[0], grid);
  if (myid==0)  printf("%4d:  df/dy     (central differences)\n", nio);

  /* xx-derivative */

  arr2D_copy_gp(A, Agp, grid);
  arr2D_deriv_xx(Agp, D, grid);
  nio = io_save_data(D[0], grid);
  if (myid==0)  printf("%4d:  d2f/dx2   (central differences)\n", nio);


  /* yy-derivative */

  arr2D_copy_gp(A, Agp, grid);
  arr2D_deriv_yy(Agp, D, grid);
  nio = io_save_data(D[0], grid);
  if (myid==0)  printf("%4d:  d2f/dy2   (central differences)\n", nio);

  /* xy-derivative */

  arr2D_copy_gp(A, Agp, grid);
  arr2D_deriv_x(Agp, D, grid);
  arr2D_copy_gp(D, Agp, grid);
  arr2D_deriv_y(Agp, D, grid);
  nio = io_save_data(D[0], grid);
  if (myid==0)  printf("%4d:  d2f/dxdy  (central differences)\n", nio);


  /* fft transform forth and back */

  arr2D_copy(A, D, grid);
  fft_test(D[0], grid);
  nio = io_save_data(D[0], grid);
  if (myid==0)  printf("%4d:  f(x,y)    FFT fransformed twice\n", nio);


  /* x-derivative */

  arr2D_copy(A, D, grid);
  fft_deriv_x(D[0], grid);
  nio = io_save_data(D[0], grid);
  if (myid==0)  printf("%4d:  df/dx     (spectral)\n", nio);

  /* y-derivative */

  arr2D_copy(A, D, grid);
  fft_deriv_y(D[0], grid);
  nio = io_save_data(D[0], grid);
  if (myid==0)  printf("%4d:  df/dy     (spectral)\n", nio);

  /* xx-derivative */

  arr2D_copy(A, D, grid);
  fft_deriv_xx(D[0], grid);
  nio = io_save_data(D[0], grid);
  if (myid==0)  printf("%4d:  d2f/dx2   (spectral)\n", nio);

  /* yy-derivative */

  arr2D_copy(A, D, grid);
  fft_deriv_yy(D[0], grid);
  nio = io_save_data(D[0], grid);
  if (myid==0)  printf("%4d:  d2f/dy2   (spectral)\n", nio);

  /* xy-derivative */

  arr2D_copy(A, D, grid);
  fft_deriv_xy(D[0], grid);
  nio = io_save_data(D[0], grid);
  if (myid==0)  printf("%4d:  d2f/dxdy  (spectral)\n", nio);

  /* laplacian */

  arr2D_copy(A, D, grid);
  fft_laplacian(D[0], grid);
  nio = io_save_data(D[0], grid);
  if (myid==0) printf("%4d:  del^2(f)  (spectral)\n", nio);


  /* spectral derivatives in 1D */

  if (np == 1) {

    memset(df0, 0, 2*N*sizeof(double));
    memset(df1, 0, 2*N*sizeof(double));
    memset(df2, 0, 2*N*sizeof(double));

    for (i=0; i<N; i++)  {
      df0[i][0] = A[N/2][i][0];
      df1[i][0] = A[N/2][i][0];
      df2[i][0] = A[N/2][i][0];
    }

    //fft_wrap_1D(df1, grid, FORWARD);
    //fft_wrap_1D(df2, grid, FORWARD);
    //fft_wrap_1D(df2, grid, BACKWARD);

    fft_deriv_x_1D (df1, grid);
    fft_deriv_xx_1D(df2, grid);

    sprintf(filename,"%s_1D.txt", ctrl.runname);
    thefile = fopen(filename, "wt");
    fprintf(thefile, "#_1.x  2.f  3.f_x  4.f_xx\n\n"); 
    for (i=0; i<N; i++) 
      fprintf(thefile, "%10.6e  %10.6e  %10.6e  %10.6e\n", 
	                i*L/N,  df0[i][0],  df1[i][0],  df2[i][0]);
    fclose(thefile);

  }


  MPI_Finalize();
  exit(0); 
}
