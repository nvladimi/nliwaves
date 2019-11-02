#include "header.h"

static  int 		 myid, np;
static  int 		 N0, Ngrids;
static  int              diagcount;

static  int              Nbins;
static  double          *fftBins;
static  double          *fftBinsTmp;
static  double          *fftBinsAll;
static  double          *weights;
static  double           L;

static  fftw_complex    *data, *fdata;

static  char             basename[80];
static  char             msg[80];

static  double           regridTh, Zup, Zdn;

extern void spectrum_average(int grid);
extern void spectrum_add();

/*----------------------------------------------------------------*/

void spectrum_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag,
		   fftw_complex *psi, fftw_complex *psihat){

  int     i, j, k;
  double  kk;

  data  = psi;
  fdata = psihat;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  strcpy(basename, ctrl->runname);

  L        =  geom->Lx;
  N0       =  geom->N0;
  Ngrids   =  geom->Ngrids;
  regridTh =  geom->regridTh;
  Zup      =  geom->regridZup;
  Zdn      =  geom->regridZdn;

  Nbins    = N0/2 * pow(2, Ngrids-1);

  //if (Nbins > N/2) Nbins = N/2;
  //if (Nbins <= 0)  Nbins = N/4;

  fftBins    = (double *)malloc(Nbins * sizeof(double) );
  fftBinsTmp = (double *)malloc(Nbins * sizeof(double) );
  fftBinsAll = (double *)malloc(Nbins * sizeof(double) );
  weights    = (double *)malloc(Nbins * sizeof(double) );


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
  
  for (k=0; k<Nbins; k++) fftBins[k] = 0;
  diagcount = 0;
}

/*----------------------------------------------------------------*/

void spectrum_compute(int grid){

  spectrum_average(grid);

  spectrum_add();

}

/*----------------------------------------------------------------*/

void spectrum_average(int grid){

  double        a, b, w;
  int           N, n, i, j, kx, ky, k;
 
  N = N0*pow(2,grid);

  n = N/np;

  w = 1./N/N;
 
  /*-- reset bins for the current snapshot --*/
  for (k=0; k<Nbins; k++) fftBinsTmp[k] = 0;
 
  /*-- sum over angle --*/
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    k = floor(sqrt(kx*kx + ky*ky)+0.5);
 
    a = fdata[N*j+i][0];
    b = fdata[N*j+i][1];

    if (k<Nbins) fftBinsTmp[k] += a*a + b*b;

  }

  /*-- normalize for current grid --*/
  for (k=0; k<Nbins; k++) fftBinsTmp[k] = fftBinsTmp[k]*w*weights[k];

}

/*----------------------------------------------------------------*/

int spectrum_regrid(int grid, int substep, double t){

  int       N, Kup, Kdn;
  int       flag;
  const int Up   = 1;
  const int Down = 2;

  /*-- no regrid for single grids, no down-grid at nonzero substeps --*/  

  if (Ngrids == 1) return(grid);
  if ((grid == Ngrids-1) && (substep != 0)) return(grid);

  /*-- each process average its spectrum --*/

  N = N0*pow(2,grid);

  Kup = Zup*N/2;
  Kdn = Zdn*N/2;
  
  spectrum_average(grid);

  /*-- master processor collects spectrum and decide to regrid --*/

  MPI_Reduce(fftBinsTmp, fftBinsAll, Nbins, 
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  flag = 0;
  if (fftBinsAll[Kup] > regridTh)             flag=Up;
  else if ((fftBinsAll[Kdn] < regridTh) && 
	   (grid > 0) && (substep == 0) )     flag=Down;

 
  MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if  ((grid == Ngrids-1) &&  (flag == Up)) return(-1);


  /*-- regrid --*/

  if (flag == Up)  {
    sprintf(msg, "#msg: regrid up       at t = %18.12e to %dx%d (%d)\n", 
	    t, N*2, N*2, grid+1);
    io_message(msg,0); 
    //sprintf(msg, "#dbg: Up:   fft[%d] = %e\n", Kup, fftBinsAll[Kup]);
    //io_message(msg,0);
    //sprintf(msg, "#dbg: Down: fft[%d] = %e\n", Kdn, fftBinsAll[Kdn]); 
    //io_message(msg,0);

    fft_regrid_up(data, grid);
    grid++;

  }

  else if (flag == Down) {
    sprintf(msg, "#msg: regrid down     at t = %18.12e to %dx%d (%d)\n", 
	    t, N/2, N/2, grid-1);
    io_message(msg,0);
    //sprintf(msg, "#dbg: Up:   fft[%d] = %e\n", Kup, fftBinsAll[Kup]);
    //io_message(msg,0);      
    //sprintf(msg, "#dbg: Down: fft[%d] = %e\n", Kdn, fftBinsAll[Kdn]); 
    //io_message(msg,0);      

    fft_regrid_dn(data, grid);
    grid--;
  }

  return(grid);
}

/*----------------------------------------------------------------*/

void spectrum_add(){

  int k;

  /*-- add spectrum from this snapshot to earlier collected --*/
  for (k=0; k<Nbins; k++) fftBins[k] += fftBinsTmp[k];


  // sprintf(msg, "#dbg: spectrum data %4d computed\n", diagcount);
  // io_message(msg,0);

  diagcount++;
}


/*----------------------------------------------------------------*/

void spectrum_output(int filecount)
{
  int    k;
  FILE      *thefile;
  char       filename[80];
  double     pi;

  pi = 2*acos(0);

  /*-- collect bin contents from all processors --*/

  MPI_Reduce(fftBins, fftBinsTmp, Nbins, 
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /*-- print out spectrum --*/
 
  if (myid == 0) {

    sprintf(filename,"%s.spc.%04d", basename, filecount);

    thefile = fopen(filename, "a");
    fprintf(thefile, "\n\n");
    fprintf(thefile, "%%\n%% Spectrum of |psi|^2:  count = %d\n", diagcount);
    fprintf(thefile, "%%\n%% 1.k  2.|psi|^2_k\n%%\n");

    for (k=0; k<Nbins; k++) fprintf(thefile, "%20.8e  %20.8e \n", 
			   k*2*pi/L, fftBinsTmp[k]/diagcount);  

    fclose(thefile); 
  }

  // sprintf(msg, "#dbg: spectrum file %04d saved\n", filecount);
  // io_message(msg,0);


  /*-- empty bins and reset count --*/

  for (k=0; k<Nbins; k++) fftBins[k] = 0;
  diagcount = 0;

}

/*----------------------------------------------------------------*/

void spectrum_output_clean(int filecount)
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
