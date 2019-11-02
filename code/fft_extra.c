#include "header.h"

static  int 		 myid, np;
static  int 		 N0;
static  double           L;
static  double           pi;

/* ---------------------------------------------------------------- */

void fft_extra_init(geom_ptr geom){

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  N0     = geom->N0;
  L      = geom->Lx;
  pi     = acos(0)*2;

}


/* ---------------------------------------------------------------- */

void fft_test(fftw_complex *data, int grid){

  fft_wrap(data, grid, FORWARD);
  fft_wrap(data, grid, BACKWARD);

}



/* ---------------------------------------------------------------- */

void fft_save_fourier(fftw_complex *data, fftw_complex *fdata, int grid){

  int N, ntot, i;

  N = N0 * pow(2, grid);
  ntot = N*N/np;

  for (i=0; i<ntot; i++) {
    fdata[i][0] = data[i][0];
    fdata[i][1] = data[i][1];
  }

}

/* ---------------------------------------------------------------- */

void fft_make_fourier(fftw_complex *data, fftw_complex *fdata, int grid){


  fft_wrap(data, grid, FORWARD);

  fft_save_fourier(data, fdata, grid);

  fft_wrap(data, grid, BACKWARD);

}

/* ---------------------------------------------------------------- */

void fft_deriv_x_1D(fftw_complex *data, int grid){

  double   c, DD, re, im;
  int      N, i, kx;

  N = N0 * pow(2,grid);

  c = 2*pi/L;

  fft_wrap_1D(data, grid, FORWARD);
 
  
  for (i=0; i<N; i++) {

    kx = i;

    if (kx > N/2) kx = kx-N;

    DD = kx*c;

    re = data[i][0];
    im = data[i][1];

    data[i][0] = - im * DD;
    data[i][1] =   re * DD;
 
  }

  fft_wrap_1D(data, grid, BACKWARD);

  
}

/* ---------------------------------------------------------------- */

void fft_deriv_xx_1D(fftw_complex *data, int grid){

  double        c, DD;
  int           N, i, kx;

  N = N0 * pow(2,grid);
  c = 4*pi*pi/(L*L);

  fft_wrap_1D(data, grid, FORWARD);
 
  for (i=0; i<N; i++) {

    kx = i;

    if (kx > N/2) kx = kx-N;

    DD = -kx*kx*c;

    data[i][0] = data[i][0] * DD;
    data[i][1] = data[i][1] * DD;

  }

  fft_wrap_1D(data, grid, BACKWARD);
}


/* ---------------------------------------------------------------- */

void fft_deriv_x(fftw_complex *data, int grid){

  double   c, DD, re, im;
  int      N, n, i, j, ij, kx;

  N = N0 * pow(2,grid);

  n = N/np;
  c = 2*pi/L;

  fft_wrap(data, grid, FORWARD);
 
  
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;

    if (kx > N/2) kx = kx-N;

    DD = kx*c;

    re = data[N*j+i][0];
    im = data[N*j+i][1];

    data[N*j+i][0] = - im * DD;
    data[N*j+i][1] =   re * DD;
 
  }

  fft_wrap(data, grid, BACKWARD);

  
}

/* ---------------------------------------------------------------- */

void fft_deriv_y(fftw_complex *data, int grid){

  double        c, DD, re, im;
  int           N, n, i, j, ky;

  N = N0 * pow(2,grid);
  n = N/np;
  c = 2*pi/L;

  fft_wrap(data, grid, FORWARD);
 
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    ky = j + myid*n;

    if (ky > N/2) ky = ky-N;

    DD = ky*c;

    re = data[N*j+i][0];
    im = data[N*j+i][1];

    data[N*j+i][0] = - im * DD;
    data[N*j+i][1] =   re * DD;
  
  }

  fft_wrap(data, grid, BACKWARD);
}

/* ---------------------------------------------------------------- */

void fft_deriv_xx(fftw_complex *data, int grid){

  double        c, DD;
  int           N, n, i, j, kx;

  N = N0 * pow(2,grid);
  n = N/np;
  c = 4*pi*pi/(L*L);

  fft_wrap(data, grid, FORWARD);
 
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;

    if (kx > N/2) kx = kx-N;

    DD = -kx*kx*c;

    data[N*j+i][0] = data[N*j+i][0] * DD;
    data[N*j+i][1] = data[N*j+i][1] * DD;

  }

  fft_wrap(data, grid, BACKWARD);
}

/* ---------------------------------------------------------------- */

void fft_deriv_yy(fftw_complex *data, int grid){

  double        c, DD;
  int           N, n, i, j, ky;

  N = N0 * pow(2,grid);
  n = N/np;
  c = 4*pi*pi/(L*L);

  fft_wrap(data, grid, FORWARD);
 
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    ky = j + myid*n;

    if (ky > N/2) ky = ky-N;

    DD = -ky*ky*c;

    data[N*j+i][0] = data[N*j+i][0] * DD;
    data[N*j+i][1] = data[N*j+i][1] * DD;
 
  }

  fft_wrap(data, grid, BACKWARD);
}


/* ---------------------------------------------------------------- */

void fft_deriv_xy(fftw_complex *data, int grid){

  double        c, DD;
  int           N, n, i, j, kx, ky;

  N = N0 * pow(2,grid);
  n = N/np;
  c = 4*pi*pi/(L*L);


  fft_wrap(data, grid, FORWARD);
 
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    DD = -kx*ky*c;

    data[N*j+i][0] = data[N*j+i][0] * DD;
    data[N*j+i][1] = data[N*j+i][1] * DD;

  }

  fft_wrap(data, grid, BACKWARD);
}

/* ---------------------------------------------------------------- */

void fft_laplacian(fftw_complex *data, int grid){

  double        c, DD;
  int           N, n, i, j, kx, ky;

  N = N0 * pow(2,grid);
  n = N/np;
  c = 4*pi*pi/(L*L);


  fft_wrap(data, grid, FORWARD);
 
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    DD = -(kx*kx + ky*ky)*c;

    data[N*j+i][0] = data[N*j+i][0] * DD;
    data[N*j+i][1] = data[N*j+i][1] * DD;

  }

  fft_wrap(data, grid, BACKWARD);
}

/* ---------------------------------------------------------------- */

void fft_davey_stewartson_phi_x(fftw_complex *data, int grid, double nu){
  // input:   data = |psi|^2
  // output:  data = phi_x

  double        c, DD;
  int           N, n, i, j, kx, ky;

  if (nu == 0) return;

  N = N0 * pow(2,grid);
  n = N/np;
  // c = 4*pi*pi/(L*L);


  fft_wrap(data, grid, FORWARD);
 
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    DD = kx*kx / (kx*kx + nu*ky*ky);

    data[N*j+i][0] = data[N*j+i][0] * DD;
    data[N*j+i][1] = data[N*j+i][1] * DD;

  }

  if (myid == 0) {
    data[0][0] = 0;
    data[0][1] = 0;
  }

  fft_wrap(data, grid, BACKWARD);
}

/* ---------------------------------------------------------------- */

void fft_blur(fftw_complex *A, fftw_complex *B, int grid, double rb)
{
  double        h0, rr, rrb, dxdx;
  double        u1, v1, u2, v2;
  int           N, n, i, j, I, J;

  N = N0 * pow(2,grid);
  n = N/np;

  /*-- create Gaussian in real space --*/

  dxdx = L*L/N/N;
  rrb = rb*rb;
  h0  = 1/(2*pi*pi*rrb); 

  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    J = j + myid*n;
    I = i;

    if (I > N/2) I = I-N;
    if (J > N/2) J = J-N;

    rr = (I*I + J*J)*dxdx;

    B[N*j+i][0] = h0*exp(-rr/rrb);
    B[N*j+i][1] = 0;    

  }


  /*-- compute convolution in Fourier space --*/

  fft_wrap(A, grid, FORWARD);
  fft_wrap(B, grid, FORWARD);

  for (i=0; i<N*n; i++) {

    u1 = A[i][0];
    v1 = A[i][1];

    u2 = B[i][0];
    v2 = B[i][1];

    B[i][0] = u1*u2 - v1*v2;
    B[i][1] = u1*v2 + v1*u2;

  }

  fft_wrap(A, grid, BACKWARD);
  fft_wrap(B, grid, BACKWARD);

}
/*----------------------------------------------------------------*/
