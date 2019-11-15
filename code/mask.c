#include "header.h"

static  int              mask;
static  int 		 myid, np;
static  int 		 N0;
static  int              grid_old;

static  double           Lx, Ly;
static  double           pi;
static  double           beta, alpha, R0;

static  double           *reMask, *imMask, **reM, **imM;
static  fftw_complex     *Psi;

static  char             msg[80];

static void mask_compute(int grid, double dt);

/* ---------------------------------------------------------------- */

void mask_init(geom_ptr geom, phys_ptr phys, fftw_complex *psi){

  int     i, j, Nx, Ny, Ngrids;

  mask = phys->mask;

  if (mask == 0) return;
  
  Psi = psi;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  N0 =  geom->N0;
  Lx =  geom->Lx;
  Ly =  geom->Ly;
  pi =  acos(0)*2;
  alpha = phys->f_alpha;
  beta  = phys->f_beta;
  R0    = phys->f_radius;

  Ngrids     =  geom->Ngrids;
  grid_old   = -1;

  /*--- create suppementary arrays ---*/

  Nx = N0*pow(2, Ngrids-1);
  Ny = Nx/np;

  reMask  =  (double *) malloc( Nx*Ny * sizeof(double));
  imMask  =  (double *) malloc( Nx*Ny * sizeof(double));
  reM =      (double **)malloc( Ny * sizeof(double *) );
  imM =      (double **)malloc( Ny * sizeof(double *) );
  
}

/* ---------------------------------------------------------------- */

void mask_compute(int grid, double dt)
{

  int i, j, nx, ny;
  double dx, x, y, rr, rx, ry, a, RR, p;
  double offset = 0;        /* 0.5 for cell-center */

  nx = N0 * pow(2,grid);
  ny = nx/np;
  dx   = Lx/nx;

  RR = R0*R0;
  p  = alpha*R0;

  for (j=0; j<ny; j++)  reM[j] = reMask + j*nx;
  for (j=0; j<ny; j++)  imM[j] = imMask + j*nx;


  for (j=0; j<ny; j++) for (i=0; i<nx; i++) {

    x = (i + offset)*dx;
    y = (myid*ny + j + offset)*dx;

    rx = fabs(x-Lx/2);
    ry = fabs(y-Ly/2);

    if (rx > 0.5*Lx) rx = Lx-rx;  
    if (ry > 0.5*Ly) ry = Ly-ry;  

    rr = rx*rx + ry*ry; 
    a = 1 - exp( -pow(rr/RR, p) );
    a = a*beta*dt;

    reM[j][i] =  cos(a);
    imM[j][i] = -sin(a);

  }

}

/* ---------------------------------------------------------------- */


void mask_apply(int grid, double dt)
{
  int N, i;
  double u,v;

  if (mask == 0) return;
  
  N = N0 * pow(2,grid);
  N = N*N/np;

  if (grid != grid_old) {
    mask_compute(grid, dt);
    grid_old = grid;
  }

  for (i=0; i<N; i++){
    u = Psi[i][0];
    v = Psi[i][1];    
    Psi[i][0] = u*reMask[i] - v*imMask[i];
    Psi[i][1] = u*imMask[i] + v*reMask[i];
  }

}


/* ---------------------------------------------------------------- */
