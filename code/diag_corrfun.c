#include "header.h"

static  int 		 myid, np;
static  int 		 N0, Ngrids;
static  int              diagcount;

static  int              Nbins;
static  double          *avgBins;
static  double          *avgBinsTmp;
static  double          *avgBinsAll;
static  double         **weights;

static  fftw_complex    *data;

static  char             basename[80];
static  char             msg[80];


extern void corrfun_average(int grid);
extern void corrfun_add();

/*----------------------------------------------------------------*/

void corrfun_init(geom_ptr geom, char *runname, fftw_complex *psi){

  int     i, j, k, m, grid, N;
  double  r;

  data  = psi;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  strcpy(basename, runname);

  N0       =  geom->N0;
  Ngrids   =  geom->Ngrids;
  Nbins    =  N0/2;

  avgBins    = (double *)malloc(Nbins * sizeof(double) );
  avgBinsTmp = (double *)malloc(Nbins * sizeof(double) );
  avgBinsAll = (double *)malloc(Nbins * sizeof(double) );

  weights    = (double **)malloc(Ngrids * sizeof(double*) );
  weights[0] = (double *)malloc(Nbins*Ngrids * sizeof(double) );

  for (grid=0; grid<Ngrids; grid++) weights[grid] = weights[0] + grid*Nbins; 

  memset(weights[0], 0, Nbins*Ngrids*sizeof(double));


  /*-- set weights --*/

  for (grid=0; grid<Ngrids; grid++) {

    N = N0*pow(2,grid);
    m = pow(2,grid);

    for (j=-N/2; j<N/2; j++) for (i=-N/2; i<N/2; i++) {
      r = sqrt(i*i + j*j);
      k = floor(r/m + 0.5);
      if (k<Nbins) weights[grid][k]++;
    }

    for (k=0; k<Nbins; k++) weights[grid][k] =  1/weights[grid][k];

  }


  /*-- empty bins and reset count --*/
  
  memset(avgBins, 0, Nbins*sizeof(double));
  diagcount = 0;
}

/*----------------------------------------------------------------*/

void corrfun_compute(int grid){

  double  u, v, w;
  int     i, N, ntot;

  N    = N0*pow(2,grid);
  ntot = N*N/np;
  w    = 1./N;

  fft_wrap(data, grid, FORWARD);

  for (i=0; i<ntot; i++) {
      u = data[i][0];
      v = data[i][1];
      data[i][0] = (u*u + v*v)*w;
      data[i][1] = 0.;
  }

  fft_wrap(data, grid, BACKWARD);

  corrfun_average(grid);

  corrfun_add();

}

/*----------------------------------------------------------------*/

void corrfun_average(int grid){

  double        r;
  int           N, n, i, j, kx, ky, k, m;
 
  N = N0*pow(2,grid);
  m = pow(2,grid);

  n = N/np;


  /*-- reset bins for the current snapshot --*/
  memset(avgBinsTmp, 0, Nbins*sizeof(double));
 
  /*-- sum over angle --*/
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    r = sqrt(kx*kx + ky*ky);
    k = floor(r/m + 0.5);
 
    if (k<Nbins) avgBinsTmp[k] += data[N*j+i][0];

  }

  /*-- normalize for current grid --*/
  for (k=0; k<Nbins; k++) avgBinsTmp[k] = avgBinsTmp[k]*weights[grid][k];

}


/*----------------------------------------------------------------*/

void corrfun_add(){

  int k;

  /*-- add corrfun from this snapshot to earlier collected --*/
  for (k=0; k<Nbins; k++) avgBins[k] += avgBinsTmp[k];


  // sprintf(msg, "#dbg: corrfun data %4d computed\n", diagcount);
  // io_message(msg,0);

  diagcount++;
}


/*----------------------------------------------------------------*/

void corrfun_output(int filecount)
{
  int    k;
  FILE      *thefile;
  char       filename[80];

  /*-- collect bin contents from all processors --*/

  MPI_Reduce(avgBins, avgBinsTmp, Nbins, 
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /*-- print out corrfun --*/
 
  if (myid == 0) {

    sprintf(filename,"%s.cf.%04d", basename, filecount);

    thefile = fopen(filename, "w");
    fprintf(thefile, "%%\n%% Auto correlation function of psi:  count = %d\n", 
	    diagcount);
    fprintf(thefile, "%%\n%% 1.k  2.corrfun \n%%\n");

    for (k=0; k<Nbins; k++) fprintf(thefile, "%6d  %20.8e \n", 
				    k, avgBinsTmp[k]/diagcount);  

    fclose(thefile); 
  }

  // sprintf(msg, "#dbg: corrfun file %04d saved\n", filecount);
  // io_message(msg,0);


  /*-- empty bins and reset count --*/

  memset(avgBins, 0, Nbins*sizeof(double));
  diagcount = 0;

}

/*----------------------------------------------------------------*/
