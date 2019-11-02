/*----------------------------------------------------------------*/
/*                                                                */
/* Small code to read in a 2D binary file (Re,Im) and write out   */
/* the text file with the ray of data starting form point i0,j0.  */
/* The ray extends to the right until the end of domain.          */
/* The domain is assumed to be square, NxN points, with size L.   */
/* Compile simply as "gcc dataray.c -o dataray.x".                */
/*                                                                */
/* usage:                                                         */
/*              dataray.x infile outfile N i0 j0 L                */
/*                                                                */
/*  where i0,j0=0,..,N-1.  (Run without arguments for reminder.)  */
/*                                                                */
/*----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

main (int argc, char **argv)
{
  int          i, i0,j0, k, N;
  FILE        *filein;
  FILE        *fileout;
  char         namein[256];
  char         nameout[256];
  char        *buff;
  double      *data;
  void        *ptr;
  double       u, v, h, L, dx;  

  if (argc>6) { 
    strcpy(namein,argv[1]); 
    strcpy(nameout,argv[2]); 
    N  = atoi(argv[3]); 
    i0 = atoi(argv[4]);
    j0 = atoi(argv[5]);
    L  = atof(argv[6]);
  }else{ 
    printf("usage: dataray infile outfile N i0 j0 L (where i0,j0=0,..,N-1) \n"); 
    exit(1);
  }

  if ((filein = fopen(namein, "rb")) == NULL){
    printf("Data file not found\n");
    exit(1);
  }


  buff = (char*) malloc( N*16 );
  fseek(filein, 0L, SEEK_SET);
  for (i=0; i<16; i++) fseek(filein, j0*N, SEEK_CUR);
  fread(buff, 1, N*16, filein);
  fclose(filein);

  ptr=buff;
  data=ptr;

  fileout = fopen (nameout, "wt");
  fprintf(fileout, "%%1.r  2.abs(Psi)   3.Re(Psi)   4.Im(Psi)\n\n");

  dx = L/N;

  for (i=i0; i<i0+N/2; i++) {
    u = data[2*i];
    v = data[2*i+1];
    h = sqrt(u*u + v*v);
    fprintf(fileout, " %16.8e  %16.8e  %16.8e  %16.8e\n", (i-i0)*dx, h, u, v);
  }

  fclose(fileout);

  free(buff);
  exit (0);
}
