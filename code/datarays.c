/*----------------------------------------------------------------*/
/*                                                                */
/* Small code to read in a 2D binary file (Re,Im) and write out   */
/* the text file with the ray of data starting form point i0,j0.  */
/* The ray extends to the right until the end of domain.          */
/* The domain is assumed to be square, NxN points, with size L.   */
/* Compile simply as "gcc datarays.c -o datarays.x".              */
/*                                                                */
/* usage:                                                         */
/*              datarays.x infile outfile N L                     */
/*                                                                */
/* (Run without arguments for reminder.)                          */
/*                                                                */
/*----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

main (int argc, char **argv)
{
  int          i, j, i0, j0, k, N;
  FILE        *filein;
  FILE        *fileout;
  char         namein[256];
  char         nameout[256];
  char        *buff_x, *buff_y;
  double      *data_x, *data_y;
  void        *ptr;
  double       ux, vx, uy, vy, L, dx;  

  if (argc>4) { 
    strcpy(namein,argv[1]); 
    strcpy(nameout,argv[2]); 
    N  = atoi(argv[3]); 
    L  = atof(argv[4]);
  }else{ 
    printf("usage: datarays infile outfile N L\n"); 
    exit(1);
  }

  if ((filein = fopen(namein, "rb")) == NULL){
    printf("Data file not found\n");
    exit(1);
  }

  i0=N/2;
  j0=N/2;

  buff_x = (char*) malloc( N*16 );
  fseek(filein, 0L, SEEK_SET);
  for (i=0; i<16; i++) fseek(filein, j0*N, SEEK_CUR);
  fread(buff_x, 1, N*16, filein);

  buff_y = (char*) malloc( N*16 );
  fseek(filein, 0L, SEEK_SET);
  fseek(filein, i0*16, SEEK_CUR);
  ptr = buff_y;
  fread(ptr, 1, 16, filein);
  for (j=1; j<N; j++) {
    ptr += 16;
    fseek(filein, 16*(N-1), SEEK_CUR);
    fread(ptr, 1, 16, filein);
  }

  fclose(filein);


  ptr=buff_x;  data_x=ptr;
  ptr=buff_y;  data_y=ptr;

  fileout = fopen (nameout, "wt");
  fprintf(fileout, "%%1.r  2.Re(Psi(x))   3.Im(Psi(x))  4.Re(Psi(y))   5.Im(Psi(y)) \n\n");

  dx = L/N;

  for (i=i0; i<i0+N/2; i++) {
    ux = data_x[2*i];
    vx = data_x[2*i+1];
    uy = data_y[2*i];
    vy = data_y[2*i+1];
    fprintf(fileout, " %16.8e  %16.8e  %16.8e  %16.8e  %16.8e\n", (i-i0)*dx, ux, vx, uy, vy);
  }

  fclose(fileout);

  //free(buff_x);
  //free(buff_y);
  exit (0);
}
