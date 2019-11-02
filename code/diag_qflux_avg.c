#include "header.h"

static  int 		 myid, np;
static  int 		 N0, Ngrids;
static  int              diagcount;

static  int              Nbins;
static  double          *qflxBins;
static  double          *qflxBinsTmp;
static  double          *qflxBinsAll;
static  double          *weights;
static  double           L, kkmax;

static  fftw_complex    *psi, *psihat, *tmphat;

static  char             basename[80];
static  char             msg[80];


extern void qflux_average(int grid);
extern void qflux_add();

/*----------------------------------------------------------------*/

void qflux_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag, phys_ptr phys,
		fftw_complex *psi_in, fftw_complex *psihat_in){

  int     i, j, k, N;
  double  kk;

  psi    = psi_in;
  psihat = psihat_in;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  strcpy(basename, ctrl->runname);

  L        =  geom->Lx;
  N0       =  geom->N0;
  Ngrids   =  geom->Ngrids;
  N        =  N0 * pow(2, Ngrids-1);

  kkmax    =  N*N;

#ifdef Q_FLUX_FILTER
  kkmax    =  phys->f_kmin * phys->f_kmin;
#endif

  Nbins    = N/2;

  //if (Nbins > N/2) Nbins = N/2;
  //if (Nbins <= 0)  Nbins = N/4;

  qflxBins    = (double *)malloc(Nbins * sizeof(double) );
  qflxBinsTmp = (double *)malloc(Nbins * sizeof(double) );
  qflxBinsAll = (double *)malloc(Nbins * sizeof(double) );
  weights     = (double *)malloc(Nbins * sizeof(double) );

  tmphat = (fftw_complex *) malloc( N*N/np * sizeof(fftw_complex));

  /*-- set weights --*/

  for (k=0; k<Nbins; k++) weights[k] =  0;

  for (j=-Nbins; j<Nbins; j++){
   for (i=-Nbins; i<Nbins; i++) {
     kk = i*i + j*j;
     k = floor(sqrt(kk) + 0.5);
     if (k<Nbins) weights[k]++;
   }
  }

  for (k=0; k<Nbins; k++) weights[k] =  1/weights[k];


  /*-- empty bins and reset count --*/
  
  for (k=0; k<Nbins; k++) qflxBins[k] = 0;
  diagcount = 0;
}

/*----------------------------------------------------------------*/

void qflux_compute(int grid){

  double        a, b, u, v, p;
  int           Nx, Ny, n, i, j;
  int           N, kx, ky, kk; 

  Nx = N0*pow(2,grid);
  Ny = N0/np;
  N  = Nx;

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


  /*-- compute Q_k = F_k * conj(psi_k) --*/

  n = 0;

  for (j=0; j<Ny; j++)  for (i=0; i<Nx; i++)   {

      kx = i;
      ky = j + myid*Ny;

      if (kx > N/2) kx = kx-N;
      if (ky > N/2) ky = ky-N;

      kk = kx*kx + ky*ky;

      if (kk < kkmax) {

	a =  tmphat[Nx*j+i][0];
	b =  tmphat[Nx*j+i][1];

	u =  psihat[Nx*j+i][0];
	v = -psihat[Nx*j+i][1];

	tmphat[Nx*j+i][0] = a*u - b*v;
	tmphat[Nx*j+i][1] = b*u + a*v;

      } else {

	tmphat[Nx*j+i][0] = 0;
	tmphat[Nx*j+i][1] = 0;
      }

    }


  /*-- transform Q_k to r-space --*/

  fft_wrap(tmphat, grid, BACKWARD);


  /*-- average over angle --*/

  qflux_average(grid);

  qflux_add();

}

/*----------------------------------------------------------------*/

void qflux_average(int grid){

  double        a, b, w;
  int           N, n, i, j, kx, ky, k;
 
  N = N0*pow(2,grid);

  n = N/np;

  w = 1./N;
   

  /*-- reset bins for the current snapshot --*/
  for (k=0; k<Nbins; k++) qflxBinsTmp[k] = 0;
 
  /*-- sum over angle --*/
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    k = floor(sqrt(kx*kx + ky*ky)+0.5);
 
    if (k<Nbins) qflxBinsTmp[k] += -tmphat[N*j+i][1];

  }

  /*-- normalize for current grid --*/
  for (k=0; k<Nbins; k++) qflxBinsTmp[k] = qflxBinsTmp[k]*w*weights[k];

}


/*----------------------------------------------------------------*/

void qflux_add(){

  int k;

  /*-- add qflux from this snapshot to earlier collected --*/
  for (k=0; k<Nbins; k++) qflxBins[k] += qflxBinsTmp[k];


  // sprintf(msg, "#dbg: qflux data %4d computed\n", diagcount);
  // io_message(msg,0);

  diagcount++;
}


/*----------------------------------------------------------------*/

void qflux_output(int filecount)
{
  int        k;
  FILE      *thefile;
  char       filename[80];
  double     pi, dx;

  pi = 2*acos(0);

  dx = L/N0;  // assuming single grid 

  /*-- collect bin contents from all processors --*/

  MPI_Reduce(qflxBins, qflxBinsTmp, Nbins, 
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /*-- print out qflux --*/
 
  if (myid == 0) {

    sprintf(filename,"%s.flx.%04d", basename, filecount);

    thefile = fopen(filename, "a");
    fprintf(thefile, "\n\n");
    fprintf(thefile, "%%\n%% Qflux of |psi|^2:  count = %d\n", diagcount);
    fprintf(thefile, "%%\n%% 1.r  2.qflux_k\n%%\n");

    for (k=0; k<Nbins; k++) fprintf(thefile, "%20.8e  %20.8e \n", 
			   k*dx, qflxBinsTmp[k]/diagcount);  

    fclose(thefile); 
  }

  // sprintf(msg, "#dbg: qflux file %04d saved\n", filecount);
  // io_message(msg,0);


  /*-- empty bins and reset count --*/

  for (k=0; k<Nbins; k++) qflxBins[k] = 0;
  diagcount = 0;

}

/*----------------------------------------------------------------*/

void qflux_output_clean(int filecount)
{
  FILE      *thefile;
  char       filename[80];
 
  if (myid == 0) {
    sprintf(filename,"%s.spc.%04d", basename, filecount);
    thefile = fopen(filename, "w");
    fclose(thefile); 
  }

}

/*----------------------------------------------------------------*/
