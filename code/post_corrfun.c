#include "header.h"

static int ntot, grid;


/* -------------------------------------------------------------- */

void post_corrfun(post_ptr post, geom_ptr geom)

{
  fftw_complex *psi;  
  fftw_complex *work;
  float        *dataout;
  int           N0, np, i, n, fn;
  int           err;
  float         u,v;

  char          filename[80];

  /*-- initialization --*/

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  N0 = geom->N0;

  n    = N0 * pow(2, geom->Ngrids-1);
  ntot = n*n/np;


  psi     = (fftw_complex *)malloc( (ntot) * sizeof(fftw_complex) );
  work    = (fftw_complex *)malloc( (ntot) * sizeof(fftw_complex) );
  dataout = (float *)malloc( (ntot) * sizeof(float) );

  post_io_init();
  fft_init(geom, work);
 
  corrfun_init(geom, post->runname, psi);


  /*-- main loop --*/

  for (fn=post->fn_beg; fn<=post->fn_end; fn+=post->fn_inc) {


    /*-- figure out grid size --*/

    sprintf(filename, "%s.psi.%04d", post->runname, fn);

    grid = io_get_grid(filename, N0);
    n    = N0 * pow(2, grid);
    ntot = n*n/np;


    /*-- read data, write text file --*/

    io_read_bin(psi, ntot*sizeof(fftw_complex), filename);

    corrfun_compute(grid);

    corrfun_output(fn);


    /*-- write binary file --*/

    if (post->bin_out) {
      for (i=0; i<ntot; i++) dataout[i] = psi[i][0];
      sprintf(filename, "%s.cfbin.%04d", post->runname, fn);
      io_write_bin(dataout, ntot*sizeof(float), filename);
    }

  }

  free(psi);
 
}

/*----------------------------------------------------------------*/

