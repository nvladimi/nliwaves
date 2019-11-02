#include "header.h"

static  int 		 myid, np;
static  int 		 N0;
static  int              diagcount;

static  int              pdfNbins;
static  double           pdfMax;
static  double          *pdfBins;
static  double          *pdfBinsAll;
static  int             *pdfBinsCnt;

static  fftw_complex    *data;

static  char             basename[80];
static  char             msg[80];

/*----------------------------------------------------------------*/

void pdf_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag, 
	      fftw_complex *psi){

  int k;

  data = psi;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  strcpy(basename, ctrl->runname);
  N0  =  geom->N0;

  pdfMax  = diag->pdfMax;
  pdfNbins = diag->pdfNbins;

  pdfBins    = (double *)malloc(3*pdfNbins * sizeof(double) );
  pdfBinsAll = (double *)malloc(3*pdfNbins * sizeof(double) );
  pdfBinsCnt = (int *)malloc(3*pdfNbins * sizeof(int) );

  for (k=0; k<3*pdfNbins; k++) pdfBins[k] = 0;
  diagcount = 0;

}

/*----------------------------------------------------------------*/

void pdf_compute(int grid)
{
  int            N, i, k, n;
  double         a, b, psi, dpsi, w;

  N = N0*pow(2,grid);
  n = N*N/np;
  w = 1./N/N;

  dpsi   =  pdfMax/pdfNbins;

  /* reset bins for the current snapshot */
  for (k=0; k<3*pdfNbins; k++) pdfBinsCnt[k] = 0;

  /* do the count */
  for (i=0; i<n; i++) {

      a   = fabs(data[i][0]);
      b   = fabs(data[i][1]);
      psi = sqrt(a*a + b*b);

      k   = floor ( psi/dpsi );
      if (k >= pdfNbins) k=pdfNbins-1;
      pdfBinsCnt[3*k]++;

      k   = floor ( a/dpsi );
      if (k >= pdfNbins) k=pdfNbins-1;
      pdfBinsCnt[3*k+1]++;

      k   = floor ( b/dpsi );
      if (k >= pdfNbins) k=pdfNbins-1;
      pdfBinsCnt[3*k+2]++;

  }

  /* add PDF from this snapshot to earlier collected PDFs */
  for (k=0; k<3*pdfNbins; k++) pdfBins[k] += w*pdfBinsCnt[k];

  // sprintf(msg, "#msg:      pdf data %4d computed\n", diagcount);
  // io_message(msg,0);

  diagcount++;

}

/*----------------------------------------------------------------*/

void pdf_output(int filecount)
{
  int    k;
  double dpsi, dpdf;
  FILE      *thefile;
  char       filename[80];

  /* collect bin contents from all processors */

  MPI_Reduce(pdfBins, pdfBinsAll, 3*pdfNbins, 
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* print out PDF */
  if (myid == 0) {

    dpsi   =  pdfMax/pdfNbins;
    dpdf   =  1./dpsi/diagcount;

    sprintf(filename,"%s.pdf.%04d", basename, filecount);

    thefile = fopen(filename, "a");
    fprintf(thefile, "\n\n");
    fprintf(thefile, "%%\n%% PDF of psi: count = %d\n",  diagcount);
    fprintf(thefile, "%%\n%% 1.bin  2.psi  3.pdf(abs)  4.pdf(Re)  5.pdf(Im)\n%%\n");

    for (k=0; k<pdfNbins; k++)
      fprintf(thefile, "%6d  %11.4e  %20.8e  %20.8e  %20.8e\n", k, (k+0.5)*dpsi, 
	      pdfBinsAll[3*k]*dpdf,  pdfBinsAll[3*k+1]*dpdf, pdfBinsAll[3*k+2]*dpdf);  

    fclose(thefile); 
  }

  // sprintf(msg, "#msg:      pdf file %04d saved\n", filecount);
  // io_message(msg,0);

  /* empty bins and reset count*/

  for (k=0; k<3*pdfNbins; k++) pdfBins[k] = 0;
  diagcount = 0;

}

/*----------------------------------------------------------------*/

void pdf_output_clean(int filecount)
{
  FILE      *thefile;
  char       filename[80];

  if (myid == 0) {
    sprintf(filename,"%s.pdf.%04d", basename, filecount);
    thefile = fopen(filename, "w");
    fclose(thefile); 
  }

}

/*----------------------------------------------------------------*/
