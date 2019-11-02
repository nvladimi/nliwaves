#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

main (int argc, char **argv)
{

  int          N, N0, stride;
  int          i, j, n, n0;
  FILE        *thefile;
  char         namein[256];
  char         nameout[256];
  double      *row;
  double      *sieve;
  long int     s;

  /*-- figure out input and output files --*/

  if (argc>3) { 
    strcpy(namein,argv[1]); 
    strcpy(nameout,argv[2]); 
    N0=atoi(argv[3]);
  }else{ 
    printf("usage: datasieve infile outfile sievesize\n"); 
    exit(1);
  }


  /*-- allocate arrays --*/

  N = gridsize(namein);

  s = 2*sizeof(double);

  row   = (double*) malloc(N*s);
  sieve = (double*) malloc(N0*N0*s);

  stride = (int)(N/N0);



  /*-- read selected rows, save selected data --*/

  thefile = fopen (namein, "rb");

  for (j=0; j<N0; j++) {

    fseek(thefile, N*s*j*stride, SEEK_SET);
    fread(row, 1, N*s, thefile); 

    for (i=0; i<N0; i++) {
      n0 = j*N0+i;
      n  = i*stride;
      sieve[2*n0] = row[2*n];
      sieve[2*n0+1] = row[2*n+1];
    }
    
  }

  fclose(thefile);


  /*-- write seived data --*/

  
  thefile = fopen (nameout, "wb");

  fwrite(sieve, 1, N0*N0*s, thefile); 

  fclose(thefile);

  free(row);
  free(sieve);

  exit (0);


}

/*-----------------------------------------------*/

int gridsize(char  *filename) {

  FILE  *thefile;
  long int      s, filesize, checksize;
  int           N;
  
  s = 2*sizeof(double);

  if ( (thefile = fopen(filename, "rb")) == NULL ) {
    printf ("\n  File \"%s\" does not exist.\n\n", filename);
    exit(1);
  }

  fseek(thefile, 0L, SEEK_END);
  filesize = ftell(thefile);
  fclose(thefile);

  N    = (int) sqrt(filesize/s);
  checksize = s*N*N;

  if (filesize != checksize){
    printf("\n  Cannot determine grid for file \"%s\"\n", filename);
    printf("    File size = %ld, N = %d\n", filesize, N);
    exit(1);
  }

  return(N);

}

/*-----------------------------------------------*/

