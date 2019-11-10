#include "header.h"

static  int 		 myid, np;
static  int 		 N0;

static  fftw_complex  *Eta, *Psi, *T1, *T2, *T3;

static void   capillary_rhs(int grid, double a, double b);    // Eta, Psi = RHS(Eta, Psi)

/* --------------------------------------------- */


int main(int argc, char **argv)
{
  fftw_complex 	**work;
 
  ctrl_str  ctrl;
  geom_str  geom;
  ic_str    ic;

  int grid = 0;

  double   L = 1;
  int      N = 128;
  int      Nx, Ny;
  int      nio;

  char   filename[80];
  FILE   *thefile;
  int     i;

  double   a=1;
  double   b=1;

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


  N0 = N;

  /* major array allocation */

  arr2D_init(&geom, 0);

  work    = arr2D_create(grid, 0);

  /*--- create suppementary arrays ---*/

  Nx = N0*pow(2, geom.Ngrids-1);
  Ny = Nx/np;

  Eta = (fftw_complex *) malloc( Nx*Ny * sizeof(fftw_complex));
  Psi = (fftw_complex *) malloc( Nx*Ny * sizeof(fftw_complex));
  T1 =  (fftw_complex *) malloc( Nx*Ny * sizeof(fftw_complex));
  T2 =  (fftw_complex *) malloc( Nx*Ny * sizeof(fftw_complex));
  T3 =  (fftw_complex *) malloc( Nx*Ny * sizeof(fftw_complex));


  /* initialization */

  fft_init(&geom, work[0]);
  io_init(&ctrl, &geom);
  ic_set(Eta, &geom, &ic);

/* ------------------------------------------------- */

  /* field variable f(x,y) */

  nio = io_save_data(Eta, grid);
  if (myid==0)  printf("%4d:  capillary Eta\n", nio);

  nio = io_save_data(Psi, grid);
  if (myid==0)  printf("%4d:  capillary Psi\n", nio);


  capillary_rhs(grid, a, b);

  nio = io_save_data(Eta, grid);
  if (myid==0)  printf("%4d:  capillary dEta/dt \n", nio);

  nio = io_save_data(Psi, grid);
  if (myid==0)  printf("%4d:  capillary dPsi/dt\n", nio);

  MPI_Finalize();
  exit(0); 
}

/* ---------------------------------------------------------------- */

void capillary_rhs(int grid, double a, double b)   // Eta, Psi = rhs(Eta, Psi) 
{

  int i, N;
 
  N = N0 * pow(2,grid);
  N = N*N/np;

  for (i=0; i<N; i++){
    T1[i][0]  = Psi[i][0];
    T1[i][1]  = 0;
    T2[i][0]  = Psi[i][0];
    T2[i][1]  = 0;
    T3[i][0]  = Eta[i][0];
    T3[i][1]  = 0;
  }

  
  fft_deriv_x(T1, grid);
  fft_deriv_y(T2, grid);
  fft_laplacian(T3, grid);

  /*--  compute dPsi = laplacian(Eta) - b ( grad(Psi) )^2 / 2 --*/ 
  /*--  compute [T1,T2] = (1 + a Eta) grad(Psi) --*/

  for (i=0; i<N; i++){

    Psi[i][0] = T3[i][0] - 0.5 * b * (T1[i][0]*T1[i][0] + T2[i][0]*T2[i][0]);
    Psi[i][1] = 0;

    T1[i][0]  = (1 + a * Eta[i][0]) * T1[i][0];
    T1[i][1]  = 0;
    
    T2[i][0]  = (1 + a * Eta[i][0]) * T2[i][0];
    T2[i][1]  = 0;
  }

  /*--  compute dEta = - div( (1 + a Eta) grad(Psi) )  --*/

  
  fft_deriv_x(T1, grid);
  fft_deriv_y(T2, grid);

  for (i=0; i<N; i++){

    Eta[i][0] = - T1[i][0] - T2[i][0];
    Eta[i][1] = 0;

  }

  
}

/* ---------------------------------------------------------------- */
