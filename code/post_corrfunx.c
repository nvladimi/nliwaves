#include "header.h"

static int            myid, np, N;
  
static double        *S2, *S4, *S6, *S8;

void corrfunx_compute(fftw_complex *psi);
void corrfunx_add();
void corrfunx_avg(int nf);
void corrfunx_init();


/* -------------------------------------------------------------- */

void post_corrfunx(post_ptr post, geom_ptr geom){

  fftw_complex  *psi;
  int           fn, nf, NN, n;
  int           err;

  char          filename[80];
  FILE         *thefile;


  N  = geom->N0;
  NN = N*N;

  post_io_init();
  corrfunx_init();

  n = NN/np;

  psi  = (fftw_complex *)malloc( NN*sizeof(fftw_complex) );


  for (fn=post->fn_beg; fn<=post->fn_end; fn+=post->fn_inc) {

    sprintf(filename, "%s.psi.%04d", post->runname, fn);

    thefile = fopen(filename, "rb");
    fread(psi, NN*sizeof(fftw_complex), 1, thefile);
    fclose(thefile); 


    corrfunx_compute(psi);


    sprintf(filename, "%s.sf2.%04d", post->runname, fn);
    io_write_bin(S2, n*sizeof(double), filename);
 
    sprintf(filename, "%s.sf4.%04d", post->runname, fn);
    io_write_bin(S4, n*sizeof(double), filename);

    sprintf(filename, "%s.sf6.%04d", post->runname, fn);
    io_write_bin(S6, n*sizeof(double), filename);

    sprintf(filename, "%s.sf8.%04d", post->runname, fn);
    io_write_bin(S8, n*sizeof(double), filename);

  }

}


/*----------------------------------------------------------------*/

void corrfunx_init(){
  int n;

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  n  = N*N/np;

  S2   = (double *)malloc( n*sizeof(double) );
  S4   = (double *)malloc( n*sizeof(double) );
  S6   = (double *)malloc( n*sizeof(double) );
  S8   = (double *)malloc( n*sizeof(double) );

}


/*----------------------------------------------------------------*/

void corrfunx_compute(fftw_complex *psi){

  int   n, Nx, Ny, NN;
  int   I, J, J0, K;
  int   i1, j1, k1, i2, j2, k2;

  double u, v, s2, s4;

  NN = N*N;
  n = NN/np; 

  Nx = N;
  Ny = N/np;

  J0 = myid*Ny;

  memset(S2, 0, n*sizeof(double));
  memset(S4, 0, n*sizeof(double));
  memset(S6, 0, n*sizeof(double));
  memset(S8, 0, n*sizeof(double));


  for (J=0; J<Ny; J++)  for (I=0; I<Nx; I++) {

    K = J*Nx + I;

    for (j1=0; j1<N; j1++) for (i1=0; i1<N; i1++) {

      
      i2 = i1 + I;
      j2 = j1 + J + J0;

      if (i2>=N) i2=i2-N;
      if (j2>=N) j2=j2-N;

      k1 = j1*N + i1;
      k2 = j2*N + i2;

      u = psi[k1][0] - psi[k2][0];
      v = psi[k1][1] - psi[k2][1];

      s2 = u*u + v*v;
      s4 = s2*s2; 

      S2[K] += s2;
      S4[K] += s4;
      S6[K] += s2*s4;
      S8[K] += s4*s4;


    }
  }

  /*-- compute averages --*/

  for (I=0; I<n; I++) {

    S2[I] = S2[I]/NN;
    S4[I] = S4[I]/NN;
    S6[I] = S6[I]/NN;
    S8[I] = S8[I]/NN;

  }


}

/*----------------------------------------------------------------*/
