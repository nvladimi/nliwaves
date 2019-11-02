#include "header.h"

static int    np, myid;
static int    filecount;
static int    sievecount;
static int    klowcount;
static int    tmrRead;
static int    tmrWrite;
static int    N0;

static int    Ns, Nk;
static int    sbuffsize;
static fftw_complex *sbuff, *kbuff;

static char   runname[80];
       char   filename[80];  
       char   msg[80];

       FILE   *thefile;

/*-----------------------------------------------------------*/


void io_init(ctrl_ptr ctrl, geom_ptr geom)
{
  int N;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  strcpy (runname,  ctrl->runname);
 
  filecount = ctrl->restart;

  tmrRead  = timer_set("io:binread");
  tmrWrite = timer_set("io:binwrite");

  N0 = geom->N0;
  Ns = geom->sieve;
  Nk = geom->klow;
  N  = N0 * pow(2, geom->Ngrids - 1);


  if (Ns) {

    sbuffsize = Ns*Ns/np;
    if (sbuffsize < Ns) sbuffsize = Ns;
    sbuffsize = sbuffsize * sizeof(fftw_complex);
    sbuff = (fftw_complex *) malloc( sbuffsize );

  }

  if (Nk) {

    kbuff = (fftw_complex *) malloc(Nk*N/np * sizeof(fftw_complex));

  }



}

/*-----------------------------------------------------------*/

int io_save_data(fftw_complex *buff, int grid)
{
  int            Nx, Ny, flag;
  long int       buffsize;
  MPI_Status    status;

  timer_on(tmrWrite);

  Nx = N0*pow(2, grid);
  Ny = Nx/np;

  buffsize = Nx*Ny * sizeof(fftw_complex);

  sprintf(filename,"%s.psi.%04d", runname, filecount);

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0) flag = 1;
  else  MPI_Recv(&flag, 1, MPI_INTEGER, myid-1, myid, MPI_COMM_WORLD, &status);


  if (myid==0) thefile = fopen (filename, "wb");
  else thefile = fopen (filename, "ab");

  fwrite (buff, buffsize, 1, thefile);

  fclose(thefile); 


  if (myid != np-1)
    MPI_Send(&flag, 1, MPI_INTEGER, myid+1, myid+1, MPI_COMM_WORLD);


  filecount++;
  sievecount = 0;
  klowcount  = 0;

  timer_off(tmrWrite);

  return(filecount-1);

}
/*-----------------------------------------------------------*/

void io_save_sieve(fftw_complex *buff, int grid)
{
  int            Nx, Ny, flag;
  int            i, j, ns, stride, idstride;

  MPI_Status    status;

  if (Ns == 0) return;

  timer_on(tmrWrite);

  /*-- figure out dimensions and strides --*/

  Nx = N0*pow(2, grid);
  Ny = Nx/np;

  ns = Ns/np;
  if (ns < 1) ns = 1;

  idstride = np/Ns;
  if (idstride < 1) idstride = 1;

  stride = Nx/Ns;

  /*-- copy selected data --*/

  for (j=0; j<ns; j++) {
    for (i=0; i<Ns; i++) {
      sbuff[j*Ns+i][0] = buff[j*stride*Nx + i*stride][0];
      sbuff[j*Ns+i][1] = buff[j*stride*Nx + i*stride][1];
    }
  }


  /*-- write selected data --*/

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0) flag = 1;
  else  MPI_Recv(&flag, 1, MPI_INTEGER, myid-1, myid, MPI_COMM_WORLD, &status);


  if (remainder(myid, idstride) == 0) {

    sprintf(filename,"%s.siv.%04d", runname, filecount);
    if (myid==0 && sievecount==0) thefile = fopen (filename, "wb");
    else thefile = fopen (filename, "ab");
    fwrite (sbuff, sbuffsize, 1, thefile);

    fclose(thefile); 

  }

  if (myid != np-1)
    MPI_Send(&flag, 1, MPI_INTEGER, myid+1, myid+1, MPI_COMM_WORLD);


  sievecount++;

  timer_off(tmrWrite);


}

/*-----------------------------------------------------------*/

void io_save_klow(fftw_complex *fA, int grid)
{

  int            Nx, Ny, flag;
  int            i, j, n, J;
  long int       NN;

  MPI_Status    status;

  if (Nk == 0) return;

  timer_on(tmrWrite);


  /*-- figure out dimensions --*/

  Nx = N0*pow(2, grid);
  Ny = Nx/np;
  NN = Nx*Nx;

  /*-- copy selected data --*/

  n = 0;  // number of points to copy;

  for (j=0; j<Ny; j++) {

    J = myid*Ny+j;       // global index

    if ((J < Nk/2) || (J >= Nx-Nk/2)) {


      // copy row
      for (i=0; i<Nk/2; i++) {
	kbuff[n][0]   = fA[j*Nx+i][0]/NN;
	kbuff[n++][1] = fA[j*Nx+i][1]/NN;
      }

      for (i=Nx-Nk/2; i<Nx; i++){
	kbuff[n][0]   = fA[j*Nx+i][0]/NN;
	kbuff[n++][1] = fA[j*Nx+i][1]/NN;
      }
    }

  }

  /*-- write selected data --*/

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0) flag = 1;
  else  MPI_Recv(&flag, 1, MPI_INTEGER, myid-1, myid, MPI_COMM_WORLD, &status);

  if (n > 0) {

    sprintf(filename,"%s.klo.%04d", runname, filecount);
    if (myid==0 && klowcount==0) thefile = fopen (filename, "wb");
    else thefile = fopen (filename, "ab");
    fwrite (kbuff, n*sizeof(fftw_complex), 1, thefile);

    fclose(thefile); 

  }

  if (myid != np-1)
    MPI_Send(&flag, 1, MPI_INTEGER, myid+1, myid+1, MPI_COMM_WORLD);


  klowcount++;

  timer_off(tmrWrite);


}

/*-----------------------------------------------------------*/

void io_save_qflux(double *fA, int nio)
{

  int            n, flag;

  MPI_Status    status;

  timer_on(tmrWrite);


  /*-- figure out dimensions --*/

  n = N0/2*N0/np;


  /*-- write data --*/

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0) flag = 1;
  else  MPI_Recv(&flag, 1, MPI_INTEGER, myid-1, myid, MPI_COMM_WORLD, &status);

  if ( (myid < np/4) || (myid >= 3*np/4) ) {

    sprintf(filename,"%s.qfx.%04d", runname, nio);
    if (myid==0) thefile = fopen (filename, "wb");
    else thefile = fopen (filename, "ab");
    fwrite (fA, n*sizeof(double), 1, thefile);

    fclose(thefile); 

  }

  if (myid != np-1)
    MPI_Send(&flag, 1, MPI_INTEGER, myid+1, myid+1, MPI_COMM_WORLD);


  timer_off(tmrWrite);


}

/*-----------------------------------------------------------*/

int io_read_data(fftw_complex *buff, int grid)
{
  int           Nx, Ny, flag, i;
  long int      buffsize;
  MPI_Status    status;

  timer_on(tmrRead);

  Nx = N0*pow(2, grid);
  Ny = Nx/np;

  buffsize = Nx*Ny * sizeof(fftw_complex);

  sprintf(filename,"%s.psi.%04d", runname, filecount);

  if (myid==0) flag = 1;
  else  MPI_Recv(&flag, 1, MPI_INTEGER, myid-1, myid, MPI_COMM_WORLD, &status);

  if ( (thefile = fopen(filename, "rb")) == NULL ) {
    printf ("\n  File \"%s\" does not exist.\n\n", filename);
    MPI_Finalize();
    exit(1);
  }else{
    fseek(thefile, 0L, SEEK_SET);
    for (i=0; i<myid; i++) fseek(thefile, buffsize, SEEK_CUR);
    // fseek(thefile, myid*buffsize, SEEK_SET);
    if (fread(buff, buffsize, 1, thefile) != 1){
      printf ("\n File \"%s\" has wrong size. \n\n", filename ); 
      MPI_Finalize();
      exit(1);
      }
    fclose(thefile); 
  }

  if (myid != np-1)
    MPI_Send(&flag, 1, MPI_INTEGER, myid+1, myid+1, MPI_COMM_WORLD);

  filecount++;
  sievecount = 0;
  klowcount  = 0;

  timer_off(tmrRead);

  return(filecount-1);

}

/*-----------------------------------------------------------*/

void io_save_tag(int nio, int n0, int grid, double t)
{

  sprintf(filename,"%s.tag", runname);

  if (myid==0) {
    thefile = fopen (filename, "wt");
    fprintf(thefile, "%d  %d  %d  # nio n0 grid \n", nio, n0, grid);  
    fclose(thefile); 
  }

  sprintf(msg, "#msg: data %04d saved at t = %18.12e\n", nio, t);
  io_message(msg,0);

}

/*-----------------------------------------------------------*/

void io_message(char *msg, int mode)
{
  sprintf(filename,"%s.out", runname);

  if ((myid == 0) || (mode == 1)) {
    thefile = fopen (filename, "a");
    fprintf(thefile,"%s", msg);
    fclose(thefile);
   }

}

/*-----------------------------------------------------------*/
