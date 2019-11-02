#include "header.h"

static int            myid, np, n, nn, M;
static fftw_complex  *psi, *psi0, *work;
static double        *psisq;  


int mode_condensate();
void psi_square_sum(fftw_complex *P, int mmax, double factor);


/* -------------------------------------------------------------- */

void post_condensate(post_ptr post, geom_ptr geom)

{
  int           N0, fn, nf;
  int           m, m0, i0;
  int           err, bin_out;
  int           grid0=0;
  long int      i;

  char          filename[80];

  /*-- initialization --*/

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  N0 = geom->N0;
  M  = post->time_slices;
  bin_out = post->bin_out;

  nn = N0*N0;
  n  = nn/np;

  psi   = (fftw_complex *)malloc( n*M*sizeof(fftw_complex) );
  psi0  = (fftw_complex *)malloc( n*sizeof(fftw_complex) );
  work  = (fftw_complex *)malloc( n*sizeof(fftw_complex) );
  psisq = (double*)malloc( n*sizeof(double) );

  post_io_init();
  fft_init(geom, work);
  fft_time_init(post);
 
  /*-- read in data --*/

  nf = 0;

  for (fn=post->fn_beg; fn<=post->fn_end; fn+=post->fn_inc) {

    sprintf(filename, "%s.psi.%04d", post->runname, fn);

    io_read_bin(&psi[nf*n], n*sizeof(fftw_complex), filename);
 
    nf++;

  }

  /*-- make Fourie transform in space  --*/

  for (m=0; m<M; m++) fft_wrap(&psi[m*n], grid0, FORWARD);


  /*-- write to file averaged spectrum --*/
  
  psi_square_sum(psi, M, 1./(N0*N0*M));

  sprintf(filename, "%s.nk", post->runname);
  io_write_bin(psisq, n*sizeof(double), filename);


  /*-- write to file instant spectra --*/

  if (bin_out) {
    nf = 0;
    for (fn=post->fn_beg; fn<=post->fn_end; fn+=post->fn_inc) {
      sprintf(filename, "%s.nk.%04d", post->runname, fn);
      psi_square_sum(&psi[nf*n], 1, 1./(N0*N0));
      io_write_bin(psisq, n*sizeof(double), filename);
      nf++;
    }
  }


  /*-- separate condensate and fluctuations --*/

  fft_time_wrap(psi, FORWARD);

  sprintf(filename, "%s.wk0", post->runname);
  m0 = mode_condensate(filename, post->dtime);
  i0 = n*m0;

  for (i=0; i<n; i++) {
    psi0[i][0] = psi[i0+i][0]/sqrt(M);
    psi0[i][1] = psi[i0+i][1]/sqrt(M);
    psi[i0+i][0] = 0;
    psi[i0+i][1] = 0;
  }

  fft_time_wrap(psi, BACKWARD);


  /*-- write to file n_k for condensate and fluctuations --*/

  psi_square_sum(psi, M, 1./(N0*N0*M));
  sprintf(filename, "%s_flc.nk", post->runname);
  io_write_bin(psisq, n*sizeof(double), filename);

  psi_square_sum(psi0, 1, 1./(N0*N0));
  sprintf(filename, "%s_cnd.nk", post->runname);
  io_write_bin(psisq, n*sizeof(double), filename);


  /*-- write to file instant n_k for fluctuations --*/

  if (bin_out) {
    nf = 0;
    for (fn=post->fn_beg; fn<=post->fn_end; fn+=post->fn_inc) {
      sprintf(filename, "%s_flc.nk.%04d", post->runname, fn);
      psi_square_sum(&psi[nf*n], 1, 1./(N0*N0));
      io_write_bin(psisq, n*sizeof(double), filename);
      nf++;
    }
  }

  /*-- transform back to r-space --*/

  for (m=0; m<M; m++) fft_wrap(&psi[m*n], grid0, BACKWARD);

  fft_wrap(psi0, grid0, BACKWARD);


  /*-- write to file condensate(r) and all fluctuations(r) --*/
 
  nf = 0;
  for (fn=post->fn_beg; fn<=post->fn_end; fn+=post->fn_inc) {
    sprintf(filename, "%s_flc.psi.%04d", post->runname, fn);
    io_write_bin(&psi[nf*n], n*sizeof(fftw_complex), filename);
    nf++;
  }

  sprintf(filename, "%s_cnd.psi", post->runname);
  io_write_bin(psi0, n*sizeof(fftw_complex), filename);

  free(psi);
 
}

/*----------------------------------------------------------------*/

int mode_condensate(char *filename, double dt){

  double        u, v, p, p0, w, pi, a;
  int           m, m0;
  FILE         *thefile;

  m0 = 0;
  p0 = 0;

  pi =  acos(0)*2;

  a = 1./M/nn;

  /*-- master processor finds dominant mode --*/

  if (myid == 0) {

    thefile = fopen(filename, "wt");

    fprintf(thefile, "%% 1.m  2.omega  3.n_k  4.Re(psi)  5.Im(psi)\n");

    for (m=0; m<M; m++) {

       u = psi[m*n][0];
       v = psi[m*n][1];
       p = u*u + v*v;
       if (p>p0) {p0=p; m0=m;}

       if (m<M/2) w = 2*pi/(M*dt)*m;
       else       w = 2*pi/(M*dt)*(m-M);

       fprintf(thefile, "%6d %12.6e %12.6e %12.6e %12.6e\n", m, w, p*a, u*a, v*a);
    }

    fclose(thefile); 

  }

  MPI_Bcast( &m0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return(m0);

}


/*----------------------------------------------------------------*/

void psi_square_sum(fftw_complex *P, int mmax, double factor) {

  int    m, i;
  double u, v;

  memset(psisq, 0, n*sizeof(double));

  for (m=0; m<mmax; m++) for (i=0; i<n; i++) {
    u = P[m*n + i][0];
    v = P[m*n + i][1];
    psisq[i] += u*u + v*v;
  }

  for (i=0; i<n; i++) psisq[i] = psisq[i]*factor;

}

/*----------------------------------------------------------------*/
