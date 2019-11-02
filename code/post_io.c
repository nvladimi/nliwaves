#include "header.h"

static int    np, myid;
static int    tmrRead;
static int    tmrWrite;

/*-----------------------------------------------------------*/


void post_io_init()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  tmrRead  = timer_set("io:binread");
  tmrWrite = timer_set("io:binwrite");

}


/*-----------------------------------------------------------*/

int io_get_grid(char  *filename, int N0) {

  FILE         *thefile;
  long int      filesize, checksize;
  int           N, grid;  


  if (myid == 0) {

    if ( (thefile = fopen(filename, "rb")) == NULL ) {
      printf ("\n  File \"%s\" does not exist.\n\n", filename);
      grid = -1;
    }

    else {

      fseek(thefile, 0L, SEEK_END);
      filesize = ftell(thefile);
      fclose(thefile);

      N    = (int) sqrt(filesize/16L);
      grid = (int) log2(1.*N/N0);
      N    = N0*pow(2,grid);
      checksize = 16L*N*N;

      if (filesize != checksize){
        printf("\n  Cannot determine grid for file \"%s\"\n", filename);
	printf("    File size = %ld, N = %d, N0 = %d\n", filesize, N, N0);
        grid = -1;
      }
    }

  }

  MPI_Bcast(&grid, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

  if (grid < 0) {
    MPI_Finalize();
    exit(1);
  }
    
  return(grid);

}

/*-----------------------------------------------------------*/

void io_read_bin(void *buff, int buffsize, char *filename)
{
  FILE         *thefile;
  int           flag;
  MPI_Status    status;

  timer_on(tmrRead);

  if ( (thefile = fopen(filename, "rb")) == NULL ) {
    printf ("\n  File \"%s\" does not exist.\n\n", filename);
    MPI_Finalize();
    exit(1);
  }else{
    fseek(thefile, myid*buffsize, SEEK_SET);
    if (fread(buff, buffsize, 1, thefile) != 1){
      printf ("\n File \"%s\" has wrong size. \n\n", filename ); 
      MPI_Finalize();
      exit(1);
      }
    fclose(thefile); 
  }

  timer_off(tmrRead);

}

/*-----------------------------------------------------------*/

void io_write_bin(void *buff, int buffsize, char *filename)
{
  FILE         *thefile;
  int           flag;
  MPI_Status    status;

  timer_on(tmrWrite);

  if (myid==0) flag = 1;
  else  MPI_Recv(&flag, 1, MPI_INTEGER, myid-1, myid, MPI_COMM_WORLD, &status);


  if (myid==0) thefile = fopen (filename, "wb");
  else thefile = fopen (filename, "ab");

  fwrite (buff, buffsize, 1, thefile);

  fclose(thefile); 


  if (myid != np-1)
    MPI_Send(&flag, 1, MPI_INTEGER, myid+1, myid+1, MPI_COMM_WORLD);


  timer_off(tmrWrite);

}


/*-----------------------------------------------------------*/

void io_read_float_gp(float *A, float *buff, 
		      int nx, int ny, int gp, char *filename)
{
  FILE         *thefile;
  int           i, j, j1, j2, jj;
  int           ch1_beg, ch1_end, ch1_size;
  int           ch2_beg, ch2_end, ch2_size;
  int           ch3_beg, ch3_end, ch3_size;
  int           Nx, Ny, linesize;
  float        *ptr;    // check it - used to be void pointer 

  timer_on(tmrRead);

  Nx = nx;
  Ny = ny*np;

  linesize = Nx*sizeof(float); 

  /*-- find beginning and size for each chunk of data --*/

  j1 = myid*ny     - gp;
  j2 = (myid+1)*ny + gp;

  if (j1>=0) {
    ch1_beg = 0;
    ch1_end = 0;
    ch2_beg = j1*linesize;
  } else {
    ch1_beg = (j1+Ny)*linesize;
    ch1_end = Ny*(linesize);
    ch2_beg = 0;
  }

  if (j2<Ny) {
    ch2_end = j2*linesize;
    ch3_beg = 0;
    ch3_end = 0;
  } else {
    ch2_end = Ny*linesize;
    ch3_beg = 0;
    ch3_end = (j2-Ny)*linesize;
  }

  ch1_size = ch1_end - ch1_beg;
  ch2_size = ch2_end - ch2_beg;
  ch3_size = ch3_end - ch3_beg;

  /*
  printf("%d:  gp = %d,  j1 = %d,  j2 = %d\n", myid, gp, j1,j2);
  printf("%d:  chunk1:  beg= %4d,  end= %4d,  size= %4d\n", myid,
	 ch1_beg/linesize,  ch1_end/linesize, ch1_size/linesize);
  printf("%d:  chunk2:  beg= %4d,  end= %4d,  size= %4d\n", myid,
	 ch2_beg/linesize,  ch2_end/linesize, ch2_size/linesize);
  printf("%d:  chunk3:  beg= %4d,  end= %4d,  size= %4d\n", myid,
	 ch3_beg/linesize,  ch3_end/linesize, ch3_size/linesize);
  */


  /*-- read data from file to temporarily buffer --*/

  if ( (thefile = fopen(filename, "rb")) == NULL ) {
    printf ("\n  File \"%s\" does not exist.\n\n", filename);
    MPI_Finalize();
    exit(1);
  }

  ptr = buff;

  fseek(thefile, ch1_beg,     SEEK_SET);
  fread(ptr,     ch1_size, 1, thefile);
  ptr = ptr +    ch1_size;

  fseek(thefile, ch2_beg,     SEEK_SET);
  fread(ptr,     ch2_size, 1, thefile);
  ptr = ptr +    ch2_size;

  fseek(thefile, ch3_beg,     SEEK_SET);
  fread(ptr ,    ch3_size, 1, thefile);
  ptr = ptr +    ch3_size;

  fclose(thefile); 
 
  //printf("done reading from file\n");
  

  /*-- fill ghost points in x-direction --*/

  for (j=0; j<ny+2*gp; j++){

    jj = j*(Nx+2*gp);

    for (i=0; i<Nx; i++)  A[jj+gp+i]    = buff[j*Nx+i];
    for (i=0; i<gp; i++)  A[jj+i]       = A[jj+Nx+i];
    for (i=0; i<gp; i++)  A[jj+gp+Nx+i] = A[jj+gp+i];

  }


  timer_off(tmrRead);

}


/*-----------------------------------------------------------*/
