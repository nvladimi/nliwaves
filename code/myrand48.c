#include "header.h"

static int    np, myid;

static char   runname[80];
       char   filename[80];  
       char   msg[80];

       FILE   *thefile;

static  unsigned short currentseed[3];


/*-----------------------------------------------------------*/

double myrand48()
{

     return(erand48(currentseed));

}
     

/*-----------------------------------------------------------*/

void myrand48_save(int filecount)
{
  int            Nx, Ny, flag;
  long int       buffsize;
  MPI_Status    status;

  
  buffsize = 3 * sizeof(unsigned short);

  sprintf(filename,"%s.seed.%04d", runname, filecount);

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0) flag = 1;
  else  MPI_Recv(&flag, 1, MPI_INTEGER, myid-1, myid, MPI_COMM_WORLD, &status);


  if (myid==0) thefile = fopen (filename, "wb");
  else thefile = fopen (filename, "ab");

  fwrite (currentseed, buffsize, 1, thefile);

  fclose(thefile); 


  if (myid != np-1)
    MPI_Send(&flag, 1, MPI_INTEGER, myid+1, myid+1, MPI_COMM_WORLD);

}


/*-----------------------------------------------------------*/

void myrand48_read(ctrl_ptr ctrl, int filecount)
{
  int           Nx, Ny, flag, i;
  long int      buffsize;
  MPI_Status    status;


  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  strcpy (runname,  ctrl->runname);
   
  buffsize = 3 * sizeof(unsigned short);

  sprintf(filename,"%s.seed.%04d", runname, filecount);

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
    if (fread(currentseed, buffsize, 1, thefile) != 1){
      printf ("\n File \"%s\" has wrong size. \n\n", filename ); 
      MPI_Finalize();
      exit(1);
      }
    fclose(thefile); 
  }

  if (myid != np-1)
    MPI_Send(&flag, 1, MPI_INTEGER, myid+1, myid+1, MPI_COMM_WORLD);

}

/*-----------------------------------------------------------*/
