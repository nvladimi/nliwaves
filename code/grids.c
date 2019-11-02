#include "header.h"

static  fftw_complex    *work;

static  int  N0;
static  int  myid, np;
static  int  id_fine1, id_fine2, id_coarse;   /* for regriding */

extern void datarow2work(fftw_complex *data, int Nx, int jd, int jw);
extern void workrow2data(fftw_complex *data, int Nx, int jd, int jw);



/* ---------------------------------------------------------------- */

void grids_init(geom_ptr geom, fftw_complex *work_external){

  int  N, Ngrids;

  /*-- get geometry information --*/ 

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  work = work_external;

  Ngrids = geom->Ngrids;
  N0     = geom->N0;
  N      = N0 * pow(2, Ngrids-1);


  /*-- find process IDs for regridding --*/


  if (myid < np/2) {

    id_fine1  = 2*myid;
    id_fine2  = 2*myid + 1;
    id_coarse = floor(myid/2);

  } else {

    id_fine1  = 2*myid - np;
    id_fine2  = 2*myid - np + 1;
    id_coarse = myid + floor((np - myid)/2);  // myid/2 + np/2

  }

  //printf("%d:  fine1=%d, fine2=%d, coarse=%d\n" 
  //	 ,myid, id_fine1, id_fine2, id_coarse);
  //MPI_Finalize();  exit(1);

}

/* ---------------------------------------------------------------- */


void fft_regrid_up(fftw_complex *data, int grid){

  MPI_Status    status;
  MPI_Request   request;

  fftw_complex  *work1, *work2;

  int           Nx, Ny, ntot, i, j;


  fft_wrap(data, grid, FORWARD);

  Nx = N0 * pow(2,grid);
  Ny = Nx/np;
  ntot = Nx*Ny;

  if (np <= 2)

    /*-- simply copy data to work --*/

    for (i=0; i<ntot; i++) {
      work[i][0] = data[i][0];
      work[i][1] = data[i][1];
    }      

  else {

     memset(work, 0, ntot*4*sizeof(fftw_complex));

    /*-- even processes send data to first half of work --*/

    work1 = work;

    if ((myid < np/4) || (myid >= 3*np/4)) {
      MPI_Irecv(work1, 2*ntot, MPI_DOUBLE, id_fine1,  11, MPI_COMM_WORLD, &request);
      //printf("%d:  waiting for data1 from  %d\n", myid, id_fine1);
    }

    if ( myid%2 == 0) {
      MPI_Send(   data,  2*ntot, MPI_DOUBLE, id_coarse, 11, MPI_COMM_WORLD);
      //printf("%d:  sending data1 to        %d\n", myid, id_coarse);
    }

    if ((myid < np/4) || (myid >= 3*np/4)) {
      MPI_Wait(&request, &status);
      //printf("%d:  data1 recieved\n", myid);
    }


    /*-- odd processes send data to second half of work --*/

    work2 = &work[ntot];

    if ((myid < np/4) || (myid >= 3*np/4)) {
      MPI_Irecv(work2, 2*ntot, MPI_DOUBLE, id_fine2,  12, MPI_COMM_WORLD, &request);
      //printf("%d:  waiting for data2 from  %d\n", myid, id_fine2);
    }

    if ( myid%2 == 1) {
      MPI_Send(   data,  2*ntot, MPI_DOUBLE, id_coarse, 12, MPI_COMM_WORLD);
      //printf("%d:  sending data2 to        %d\n", myid, id_coarse);
    }

    if ((myid < np/4) || (myid >= 3*np/4)) {
      MPI_Wait(&request, &status);
      //printf("%d:  data2 recieved\n", myid);
    }
  }


  /*-- copy work to data --*/


  memset(data, 0, ntot*4*sizeof(fftw_complex));

  switch (np) {

  case 1:

    for (j=0;      j<Ny/2; j++)   workrow2data(data, Nx, j, j);
    for (j=Ny/2;   j<Ny;   j++)   workrow2data(data, Nx, j+Ny, j);
    break;

  case 2:

    if (myid == 0)
      for (j=0; j<Ny; j++)       workrow2data(data, Nx, j, j);
    else {
      for (j=0; j<Ny; j++)       workrow2data(data, Nx, j+Ny, j);
    }

    break;

  default:

    if ((myid < np/4) || (myid >= 3*np/4))
      for (j=0; j<Ny*2; j++)     workrow2data(data, Nx, j, j);
  
  }

  //io_save_data(work, grid);
  //io_save_data(data, grid+1);
  //MPI_Finalize(); exit(1);

  fft_wrap(data, grid+1, BACKWARD);

}

/*----------------------------------------------------------------*/

void fft_regrid_dn(fftw_complex *data, int grid){

  MPI_Status    status;
  MPI_Request   request;

  fftw_complex  *work1, *work2;

  int           Nx, Ny, ntot, i, j;


  fft_wrap(data, grid, FORWARD);

  Nx = N0 * pow(2,grid-1);
  Ny = Nx/np;
  ntot = Nx*Ny;


  /*-- copy data to work --*/

  memset(work, 0, ntot*4*sizeof(fftw_complex));

 
  switch (np) {

  case 1:

    for (j=0;      j<Ny/2; j++)        datarow2work(data, Nx,j,j);
    for (j=Ny/2;   j<Ny;   j++)        datarow2work(data, Nx,j+Ny,j);
    break;

  case 2:

    if (myid == 0)
      for (j=0; j<Ny; j++)             datarow2work(data, Nx,j,j);
    else
      for (j=0; j<Ny; j++)             datarow2work(data, Nx,j+Ny,j);
    break;

  default:

    if ((myid < np/4) || (myid >= 3*np/4))
      for (j=0; j<Ny*2; j++)           datarow2work(data, Nx,j,j);
  
  }



  if (np <= 2)

    /*-- simply copy work to data --*/

    for (i=0; i<ntot; i++) {
      data[i][0] = work[i][0];
      data[i][1] = work[i][1];
    }

  else {

    memset(data, 0, ntot*4*sizeof(fftw_complex));


    /*-- first half of work is sent to data on even processes --*/

    work1 = work;

    if ( myid%2 == 0) {
      MPI_Irecv(data, 2*ntot, MPI_DOUBLE, id_coarse,  13, MPI_COMM_WORLD, &request);
      //printf("%d:  waiting for data from   %d\n", myid, id_fine1);
    }

    if ((myid < np/4) || (myid >= 3*np/4)) {
      MPI_Send(work1,  2*ntot, MPI_DOUBLE, id_fine1, 13, MPI_COMM_WORLD);
      //printf("%d:  sending work1 to        %d\n", myid, id_fine1);
    }

    if ( myid%2 == 0) {
      MPI_Wait(&request, &status);
      //printf("%d:  work1 recieved\n", myid);
    }


    /*-- second half of work is sent to data on odd processes --*/

    work2 = &work[ntot];

    if ( myid%2 == 1) {
      MPI_Irecv(data, 2*ntot, MPI_DOUBLE, id_coarse,  14, MPI_COMM_WORLD, &request);
      //printf("%d:  waiting for data from   %d\n", myid, id_fine1);
    }

    if ((myid < np/4) || (myid >= 3*np/4)) {
      MPI_Send(work2, 2*ntot, MPI_DOUBLE, id_fine2, 14, MPI_COMM_WORLD);
      //printf("%d:  sending work2 to        %d\n", myid, id_fine2);
    }

    if ( myid%2 == 1) {
      MPI_Wait(&request, &status);
      //printf("%d:  work2 recieved\n", myid);
    }
  }

  fft_wrap(data, grid-1, BACKWARD);

}

/*----------------------------------------------------------------*/

void workrow2data(fftw_complex *data, int Nx, int jd, int jw) {

  int i, ij_data, ij_work;

  for (i=0;    i<Nx/2; i++) {
    ij_data = jd*Nx*2 + i;
    ij_work = jw*Nx+i;
    data[ij_data][0] = 2*work[ij_work][0];
    data[ij_data][1] = 2*work[ij_work][1];
  }
  for (i=Nx/2; i<Nx;   i++) {
    ij_data = jd*Nx*2 + Nx + i;
    ij_work = jw*Nx + i;
    data[ij_data][0] = 2*work[ij_work][0];
    data[ij_data][1] = 2*work[ij_work][1];
  }
}


void datarow2work(fftw_complex *data, int Nx, int jd, int jw) {

  int i, ij_data, ij_work;

  for (i=0;   i<Nx/2; i++) {
    ij_data = jd*Nx*2 + i;
    ij_work = jw*Nx+i;
    work[ij_work][0] = data[ij_data][0] / 2;
    work[ij_work][1] = data[ij_data][1] / 2;
  }
  for (i=Nx/2; i<Nx;   i++) {
    ij_data = jd*Nx*2 + Nx + i;
    ij_work = jw*Nx + i;
    work[ij_work][0] = data[ij_data][0] / 2;
    work[ij_work][1] = data[ij_data][1] / 2;
  }

}

/*----------------------------------------------------------------*/
