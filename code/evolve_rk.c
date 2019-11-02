#include "header.h"

static  int 		 myid, np;
static  int 		 N0;
static  double           L;
static  double           pi;

static  double           dealiasZ;
static  int              filtersize;
static  double           *filter;


static  fftw_complex  *Psi, *PsiHat, *S, *S0, *T1, *T2;

static  double           coefRho, coefNu;

static  int              tmrRK;
static  char             msg[80];


static void   compute_rhs(int grid);                  // S = RHS(S)
static void   dealias(int grid);

static void   rk_zero_out(int N);
static void   rk_add_slope(int N, int w);
static void   rk_part_step(int N, double dt);
static void   rk_full_step(int N, double dt);

/* ---------------------------------------------------------------- */

void
evolve_init(geom_ptr geom, phys_ptr phys, 
	       fftw_complex *psi, fftw_complex *psihat){

  int     i, j, Nx, Ny, filtersize, Ngrids;
  double  x;

  Psi = psi;
  PsiHat = psihat;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  N0 =  geom->N0;
  L  =  geom->Lx;
  pi =  acos(0)*2;

  coefRho = phys->coefRho;
  coefNu  = phys->coefNu;

  dealiasZ   =  geom->dealiasZ;
  filtersize =  geom->dealiasF;
  Ngrids     =  geom->Ngrids;

  /*--- set up dealising filter ---*/

  filter = (double *)malloc( filtersize*sizeof(double) );

  for (i=0; i<filtersize; i++){
    x = 1 - 1.*i/filtersize;
    filter[i] = x*x*(3-2*x);
  }

  /*--- create suppementary arrays ---*/

  Nx = N0*pow(2, Ngrids-1);
  Ny = Nx/np;

  S  =  (fftw_complex *) malloc( Nx*Ny * sizeof(fftw_complex));
  S0 =  (fftw_complex *) malloc( Nx*Ny * sizeof(fftw_complex));
  T1 =  (fftw_complex *) malloc( Nx*Ny * sizeof(fftw_complex));
  T2 =  (fftw_complex *) malloc( Nx*Ny * sizeof(fftw_complex));

 
  /*--- timers ---*/

  tmrRK    = timer_set("RK4");
 
}

/* ---------------------------------------------------------------- */


void evolve_one_step(int grid, double dt)   // Uo := Uo + dU
{
  int N;

  timer_on(tmrRK);

  N = N0 * pow(2,grid);
  N = N*N/np;

  rk_zero_out(N);

  /*-- stage 1 --*/

  rk_part_step(N, 0);
  compute_rhs(grid);
  rk_add_slope(N, 1);


  /*-- stage 2 --*/

  rk_part_step(N, dt/2);
  compute_rhs(grid);
  rk_add_slope(N, 2);


  /*-- stage 3 --*/

  rk_part_step(N, dt/2);
  compute_rhs(grid);
  rk_add_slope(N, 2);


  /*-- stage 4 --*/

  rk_part_step(N, dt);
  compute_rhs(grid);
  rk_add_slope(N, 1);


 /*-- update solution --*/

  rk_full_step(N, dt/6);

  dealias(grid);

  timer_off(tmrRK);

}

/*--------------------------------------------------------------*/

void rk_zero_out(int N)              // S=0,  S0=0;
{
  memset(S,  0, N*sizeof(fftw_complex));
  memset(S0, 0, N*sizeof(fftw_complex));
}

void rk_part_step(int N, double dt)  // U = Uo + S*h,   S: slope -> solution
{
  int i;
  for (i=0; i<N; i++) { 
    S[i][0] = Psi[i][0] + S[i][0]*dt; 
    S[i][1] = Psi[i][1] + S[i][1]*dt; 
  }
}

//void  compute_rhs(int N)          // S: solution -> slope 

void rk_add_slope(int N, int w)     // S0 = S0 + S*w
{
  int i;
  for (i=0; i<N; i++) { 
    S0[i][0] += S[i][0] * w; 
    S0[i][1] += S[i][1] * w; 
  }
}

void rk_full_step(int N, double dt)  // Uo = Uo + So*h, 
{
  int i;
  for (i=0; i<N; i++) { 
    Psi[i][0] += S0[i][0]*dt; 
    Psi[i][1] += S0[i][1]*dt; 
  }
}

/* ---------------------------------------------------------------- */

void dealias(int grid){

  int            N, n, i, j, kx, ky, dkx, dky;
  int            kmin, kmax;

  fftw_complex   psi;

  N = N0 * pow(2,grid);

  if ((dealiasZ <= 0) || (dealiasZ >= 1))  kmin = N/2;
  else                                     kmin = N/2*dealiasZ;
  kmax = kmin + filtersize;
  if (kmax > N/2) kmax = N/2; 


  n = N/np;
 
  fft_wrap(Psi, grid, FORWARD);
 
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    if ((abs(kx) < kmax) && (abs(ky) < kmax)) {

      psi[0] = Psi[N*j+i][0];
      psi[1] = Psi[N*j+i][1];
 
      /*--- smoothing dealiasing ---*/

      dkx = abs(kx) - kmin;
      dky = abs(ky) - kmin;

      if (dkx >= 0) {
	psi[0] = psi[0] * filter[dkx];
	psi[1] = psi[1] * filter[dkx];
      }

      if (dky >= 0) {
	psi[0] = psi[0] * filter[dky];
	psi[1] = psi[1] * filter[dky];
      }

    } else {
      psi[0] = 0;
      psi[1] = 0;
    }

    Psi[N*j+i][0] = psi[0];
    Psi[N*j+i][1] = psi[1];

  }

  /*-- save a copy of psihat to work array --*/

  fft_save_fourier(Psi, PsiHat, grid);

  fft_wrap(Psi, grid, BACKWARD);


}
/* ---------------------------------------------------------------- */

void evolve_biglin(int grid,  double dt){
}

/* ---------------------------------------------------------------- */


void compute_rhs(int grid)   // S = rhs(S) 
{

  int i, N;
  double u, v, pp;

  N = N0 * pow(2,grid);
  N = N*N/np;

  for (i=0; i<N; i++){
    u = S[i][0];
    v = S[i][1];
    T1[i][0]  = u;
    T1[i][1]  = v;
    T2[i][0]  = u*u + v*v;
    T2[i][1]  = 0;
  }

  fft_laplacian(T1, grid);
  fft_davey_stewartson_phi_x(T2, grid, coefNu);

  for (i=0; i<N; i++){

    u    = S[i][0];
    v    = S[i][1];
    pp   = u*u + v*v;

    S[i][0] = - T1[i][1] - pp*v  + coefRho * (u*T2[i][1] + v*T2[i][0]);
    S[i][1] =   T1[i][0] + pp*u  - coefRho * (u*T2[i][0] - v*T2[i][1]);

  }

}

/* ---------------------------------------------------------------- */
