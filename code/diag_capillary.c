#include "header.h"

void  diag_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag, phys_ptr phys,
		fftw_complex *psi, fftw_complex *psihat){

  //pdf_init       (ctrl, geom, diag, psi);
  spectrum_init  (ctrl, geom, diag, psi, psihat);
  //qflux_init     (ctrl, geom, diag, phys, psi, psihat);

}

void diag_ic(fftw_complex *psi, fftw_complex *psihat,
	     int grid, int nio, double t) {

  fft_make_fourier(psi, psihat, grid);

  info_output(grid, t);

  if (nio == 0) {

    //pdf_compute(grid);
    //pdf_output(nio);

    spectrum_compute(grid);
    spectrum_output(nio);

    //qflux_compute(grid);
    //qflux_output(nio);

  }

  //pdf_output_clean(nio+1);
  spectrum_output_clean(nio+1);
  //qflux_output_clean(nio+1);

}


void diag_compute(fftw_complex *psi, fftw_complex *psihat, int grid){

  fft_make_fourier(psi, psihat, grid);
  //pdf_compute(grid);
  spectrum_compute(grid);
  //qflux_compute(grid);
  //io_save_sieve(psi, grid);
  io_save_klow(psihat, grid);

}


void diag_output(int nio){

  //pdf_output(nio);
  spectrum_output(nio);
  //qflux_output(nio);

}

