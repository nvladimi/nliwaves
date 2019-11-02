#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

long int filesize(char  *filename);

main (int argc, char **argv)
{
  double      *A;
  double      *Atot;
  long int     i, n, n0;
  int          fn, fnum1, fnum2, fstep, nf;
  int          first_call;

  FILE        *thefile;
  char         fbase[256];
  char         filename[256];

  /*-- figure out input and output files --*/

  if (argc == 5) { 
    strcpy(fbase,argv[1]); 
    fnum1=atoi(argv[2]);
    fnum2=atoi(argv[3]);
    fstep=atoi(argv[4]);
  }else{ 
    printf("usage: dataavg  fbase fnum1  fnum2 fstep\n"); 
    exit(1);
  }


  first_call = 1;
  nf         = 0;


  /*-- loop over files --*/

  for (fn=fnum1; fn<=fnum2; fn+=fstep) {

    sprintf(filename, "%s.%04d", fbase, fn);

    n = filesize(filename);

    printf("%s: %d\n", filename, (int) n);

    /*-- if first call, allocate arrays --*/

    if (first_call) {

      n0   = n; 
      A    = (double*) malloc(n*sizeof(double));
      Atot = (double*) malloc(n*sizeof(double));
      memset(Atot, 0, n*sizeof(double));
      first_call = 0;

    } else {

      if (n != n0) {
	printf("Filesize does not match. Exiting.\n");
        exit(0);
      }

    }

    /*-- read in data and add to storage --*/

    thefile = fopen (filename, "rb");
    fread(A, 1, n*sizeof(double), thefile);
    fclose(thefile);


    for (i=0; i<n; i++) Atot[i] += A[i];
    nf++;

  }

  /*-- average data and write to file --*/

  for (i=0; i<n; i++) Atot[i] = Atot[i]/nf;


  sprintf(filename, "%s.avg", fbase);

  thefile = fopen (filename, "wb");
  fwrite(Atot, 1, n*sizeof(double), thefile);
  fclose(thefile);

  exit(0);

}

/*-----------------------------------------------*/

long int filesize(char  *filename) {

  FILE  *thefile;
  long int      s;
  
  if ( (thefile = fopen(filename, "rb")) == NULL ) {
    printf ("\n  File \"%s\" does not exist.\n\n", filename);
    exit(1);
  }

  fseek(thefile, 0L, SEEK_END);
  s = ftell(thefile) / sizeof(double);
  fclose(thefile);

  return(s);

}

/*-----------------------------------------------*/

