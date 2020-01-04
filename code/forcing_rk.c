#include "header.h"


static  fftw_complex  *Psi, *PsiHat;

static  double           dealiasZ;
static  int              filtersize;
static  double           *filter;

static  int 		 np, myid;
static  int 		 N0;


static void dealias(int grid);

/* ---------------------------------------------------------------- */


void forcing_init(geom_ptr geom, phys_ptr phys, fftw_complex *psi, fftw_complex *psihat)
{

  int     i, j, Nx, Ny, filtersize, Ngrids;
  double  x;

  Psi = psi;
  PsiHat = psihat;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  N0 =  geom->N0;
 
  dealiasZ   =  geom->dealiasZ;
  filtersize =  geom->dealiasF;
  N0 =          geom->N0;

  /*--- set up dealising filter ---*/

  filter = (double *)malloc( filtersize*sizeof(double) );

  for (i=0; i<filtersize; i++){
    x = 1 - 1.*i/filtersize;
    filter[i] = x*x*(3-2*x);
  }

 
}


/* ---------------------------------------------------------------- */

void forcing(int grid){


  fft_wrap(Psi, grid, FORWARD);

  dealias(grid);

  fft_save_fourier(Psi, PsiHat, grid);
  fft_wrap(Psi, grid, BACKWARD);
   
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

}

/* ---------------------------------------------------------------- */
