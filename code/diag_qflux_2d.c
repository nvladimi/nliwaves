#include "header.h"

static  int 		 myid, np;
static  int 		 N0, Ngrids;
static  int              diagcount;

static  fftw_complex   *tmphat;        // scratch variable
static  double         *qflux;         // storage 

static  char             basename[80];
static  char             msg[80];


/*----------------------------------------------------------------*/

void qflux_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag){

  int     N;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  strcpy(basename, ctrl->runname);

  N0       =  geom->N0;
  Ngrids   =  geom->Ngrids;

  if (Ngrids > 1) return;

  N = N0 * pow(2, Ngrids-1);

  tmphat = (fftw_complex *) malloc( N*N/np * sizeof(fftw_complex));
  qflux  = (double *) malloc( N0/2*N0/np * sizeof(double));

  memset(qflux, 0, N0/2*N0/np*sizeof(double));
  diagcount = 0;

}

/*----------------------------------------------------------------*/

void qflux_compute(fftw_complex *psi, fftw_complex *psihat, int grid){

  double        a, b, u, v, p;
  int           Nx, Ny, n, i, j;
 

  if (Ngrids > 1) return;

  Nx = N0;  //N0*pow(2,grid);
  Ny = N0/np;

  /*-- compute F=psi*|psi|^2  in r-space --*/

  for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) {

    a = psi[Nx*j+i][0];
    b = psi[Nx*j+i][1];

    p = a*a + b*b;

    tmphat[Nx*j+i][0] = a*p;
    tmphat[Nx*j+i][1] = b*p;

  }

  /*-- transform F to k-space --*/

  fft_wrap(tmphat, grid, FORWARD);


  /*-- compute Q_k = 2*Im(F_k * conj(psi_k)) --*/

  n = 0;

  for (j=0; j<Ny; j++) {

    for (i=0; i<Nx/4; i++) {

      a =  tmphat[Nx*j+i][0];
      b =  tmphat[Nx*j+i][1];

      u =  psihat[Nx*j+i][0];
      v = -psihat[Nx*j+i][1];

      qflux[n++] += 2*(b*u + a*v);
    }


    for (i=Nx-Nx/4; i<Nx; i++) {

      a =  tmphat[Nx*j+i][0];
      b =  tmphat[Nx*j+i][1];

      u =  psihat[Nx*j+i][0];
      v = -psihat[Nx*j+i][1];

      qflux[n++] += 2*(b*u + a*v);
    }

  }

  diagcount++;
}

/*----------------------------------------------------------------*/


void qflux_save(int nio){

  int           n, i;
 
  if (Ngrids > 1) return;

  n = N0*N0/2/np;

  /*-- compute F=psi*|psi|^2  in r-space --*/

  for (i=0; i<n; i++)  qflux[i] = qflux[i]/diagcount;

  io_save_qflux(qflux, nio);

  memset(qflux, 0, N0/2*N0/np*sizeof(double));
  diagcount = 0;


}


/*----------------------------------------------------------------*/






/*----------------------------------------------------------------*/

