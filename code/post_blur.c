#include "header.h"

static int ntot, grid;

extern void compute_energy(fftw_complex *P, 
			   fftw_complex *E, fftw_complex *E1);
extern void compute_magnitude(fftw_complex *P);

/* -------------------------------------------------------------- */

void post_blur(post_ptr post, geom_ptr geom)

{

  fftw_complex 	*psi, *ener, *blur, *work;
  float         *dataout;

  int           N0, np, i, n, fn;
  int           err;
  float         rb;

  char          filename[80];

  /*-- initialization --*/

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  N0 = geom->N0;

  n    = N0 * pow(2, geom->Ngrids-1);
  ntot = n*n/np;

  psi      = (fftw_complex *)malloc( (ntot) * sizeof(fftw_complex) );
  ener     = (fftw_complex *)malloc( (ntot) * sizeof(fftw_complex) );
  blur     = (fftw_complex *)malloc( (ntot) * sizeof(fftw_complex) );
  work     = (fftw_complex *)malloc( (ntot) * sizeof(fftw_complex) );
  dataout  = (float *)malloc( (ntot) * sizeof(float) );

  fft_init(geom, work);
  post_io_init();


  /*-- main loop --*/

  for (fn=post->fn_beg; fn<=post->fn_end; fn+=post->fn_inc) {


    /*-- figure out grid size --*/

    sprintf(filename, "%s.psi.%04d", post->runname, fn);

    grid = io_get_grid(filename, N0);
 
    n    = N0 * pow(2, grid);
    ntot = n*n/np;


    /*-- read data --*/

    io_read_bin(psi, ntot*sizeof(fftw_complex), filename);

    compute_energy(psi, ener, blur);  /* use blur as a scratch array */
    compute_magnitude(psi);           /* psi is overwritten */


    /*-- write un-blurred data --*/

    for (i=0; i<ntot; i++) dataout[i] = psi[i][0];
    sprintf(filename, "%s.Nblur%04d.%04d", post->runname, 0, fn);
    io_write_bin(dataout, ntot*sizeof(float), filename);

    for (i=0; i<ntot; i++) dataout[i] = psi[i][0];
    sprintf(filename, "%s.Hblur%04d.%04d", post->runname, 0, fn);
    io_write_bin(dataout, ntot*sizeof(float), filename);


    /*-- blur with each filter --*/

    for (n = 0; n < post->nblur ; n++){

      rb = post->rblur[n];

      fft_blur(psi, blur, grid, rb);

      for (i=0; i<ntot; i++) dataout[i] = blur[i][0];
      sprintf(filename, "%s.Nblur%04d.%04d", post->runname, (int)(rb*1000), fn);
      io_write_bin(dataout, ntot*sizeof(float), filename);

      fft_blur(ener, blur, grid, rb);
      for (i=0; i<ntot; i++) dataout[i] = blur[i][0];
      sprintf(filename, "%s.Hblur%04d.%04d", post->runname, (int)(rb*1000), fn);
      io_write_bin(dataout, ntot*sizeof(float), filename);

    }

  }

  free(psi);
  free(ener);
  free(blur);
  free(dataout);

}

/*----------------------------------------------------------------*/

void compute_energy(fftw_complex *P, fftw_complex *E, fftw_complex *E1)
{
  double        u, v, ux, uy, vx, vy, Ek, Ep;
  int           i;

  for (i=0; i<ntot; i++) {
    E[i][0] = P[i][0];
    E[i][1] = P[i][1];
  }

  fft_deriv_x(E, grid);

  for (i=0; i<ntot; i++) {
    E1[i][0] = P[i][0];
    E1[i][1] = P[i][1];
  }

  fft_deriv_y(E1, grid);

  for (i=0; i<ntot; i++) {

    u  = P[i][0];
    v  = P[i][1];

    ux = E[i][0];
    vx = E[i][1];
    uy = E1[i][0];
    vy = E1[i][1];

    Ek = ux*ux + vx*vx + uy*uy + vy*vy;
    Ep = u*u + v*v;
    Ep = Ep*Ep;

    E[i][0] = Ek - Ep;
    E[i][1] = 0;

  }
}

/*----------------------------------------------------------------*/

void compute_magnitude(fftw_complex *P)
{
  double        u, v;
  int           i;

  for (i=0; i<ntot; i++) {

    u  = P[i][0];
    v  = P[i][1];

    P[i][0] = u*u + v*v;
    P[i][1] = 0;

  }
}

/*----------------------------------------------------------------*/

