#include "header.h"

void  diag_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag, phys_ptr phys, 
		fftw_complex *psi, fftw_complex *psihat) {

  spectrum_init  (ctrl, geom, diag, psi, psihat);

}


void diag_ic(fftw_complex *psi, fftw_complex *psihat,
	     int grid, int nio, double t) {

  fft_make_fourier(psi, psihat, grid);
  info_output(grid, t);
  spectrum_compute(grid);

}


void diag_compute(fftw_complex *psi, fftw_complex *psihat, int grid){
}


void diag_output(int nio){
}

