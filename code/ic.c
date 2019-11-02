#include "header.h"

#define  SpBump   1
#define  SpTail   2
#define  SpExp    3
#define  SpGauss  4

static double  dx;
static double  Lx, Ly;
static int     nx, ny, gp, grid;
static int     myid, np;

void Tophat(fftw_complex **A, double x, double y, double r, double h);
void Vortex(fftw_complex **A, double x, double y, double r, double h);
void Gauss(fftw_complex **A, double x, double y, double r, double h, double phase);
void superGauss(fftw_complex **A, double x, double y, double r, double h, double p, double phase);
void Ellipse(fftw_complex **A, double x, double y, double a, double b, double h);
void MultiGauss(fftw_complex **A, double rmax, double hmax, int nbumps, int seed);
void GaussArray(fftw_complex **A, double numParticles, double rmax, int nbumps, int seed);
void GaussHex(fftw_complex **A, double numParticles, double r, double c, int m, int seed);
void superGaussHex(fftw_complex **A, double numParticles, 
                    double d, double c, double p, int m, double chirp, double width, 
                    double phaserand, double seed);

void MakeHole(fftw_complex **A, double x, double y, double r);

void icSpectrum(fftw_complex **A, double kmin, double kmax, 
	      double numParticles, double nCondensate, int seed, int SpType);

void icGaussNoise(fftw_complex **A, double kmin, double kmax, 
		  double numParticles, double nCondensate, int seed,
                  double BumpX, double BumpY, 
		  double BumpRadius, double BumpHeight);

double HoleRadius(double r, double h, double b);

void  normalizeField(fftw_complex **A, double numParticles);

/*-----------------------------------------------------------*/

void ic_set(fftw_complex *psi, geom_ptr geom, ic_ptr ic)
{

  fftw_complex **A;
 
  int j;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  grid = geom->grid;
  nx   = (geom->N0)*pow(2, grid);
  ny   = nx/np;
  Lx   = geom->Lx;
  Ly   = geom->Ly;
  dx   = Lx/nx;


  /*-- set arrays for double indexing --*/

  gp   = 0;

  A = (fftw_complex **)malloc( (ny+2*gp) * sizeof(fftw_complex *) );

  for (j=0; j<ny+2*gp; j++)  A[j] = psi + j*(nx+2*gp);


  /*-- select and IC and fill the array --*/

  if( strcmp(ic->type,"tophat") == 0 )

    Tophat(A, ic->BumpX, ic->BumpY, ic->BumpRadius, ic->BumpHeight);

  else if( strcmp(ic->type,"gauss") == 0 )

    Gauss(A, ic->BumpX, ic->BumpY, ic->BumpRadius, ic->BumpHeight, 0);

  else if( strcmp(ic->type,"ellipse") == 0 )

    Ellipse(A, ic->BumpX, ic->BumpY, ic->BumpA, ic->BumpB, ic->BumpHeight);

  else if( strcmp(ic->type,"vortex") == 0 )

    Vortex(A, ic->BumpX, ic->BumpY, ic->BumpRadius, ic->BumpHeight);

  else if( strcmp(ic->type,"multigauss") == 0 )

    MultiGauss(A, ic->BumpRadius, ic->BumpHeight, ic->NumberOfBumps, ic->Seed);

  else if( strcmp(ic->type,"gaussarray") == 0 )

    GaussArray(A, ic->numParticles, ic->BumpRadius, ic->NumberOfBumps, ic->Seed);

  else if( strcmp(ic->type,"gausshex") == 0 )

    GaussHex(A, ic->numParticles, ic->BumpRadius, ic->BumpA, ic->NumberOfBumps, ic->Seed);

  else if( strcmp(ic->type,"sghex") == 0 )

    superGaussHex(A, ic->numParticles, ic->BumpRadius, ic->BumpA, ic->BumpB, ic->NumberOfBumps, 
                     ic->Chirp, ic->BumpX, ic->BumpY, ic->Seed);


  else if( strcmp(ic->type,"sp_bump") == 0 )

    icSpectrum(A, ic->kmin, ic->kmax, ic->numParticles, 
	       ic->nCondensate, ic->Seed, SpBump);

  else if( strcmp(ic->type,"sp_tail") == 0 )

    icSpectrum(A, ic->kmin, ic->kmax, ic->numParticles, 
	       ic->nCondensate, ic->Seed, SpTail);

  else if( strcmp(ic->type,"sp_exp") == 0 )

    icSpectrum(A, ic->kmin, ic->kmax, ic->numParticles,
	       ic->nCondensate, ic->Seed, SpTail);

  else if( strcmp(ic->type,"sp_gauss") == 0 )

    icSpectrum(A, ic->kmin, ic->kmax, ic->numParticles,
	       ic->nCondensate, ic->Seed, SpTail);

  else if( strcmp(ic->type,"gauss_noise") == 0 )

    icGaussNoise(A, ic->kmin, ic->kmax, ic->numParticles,  ic->nCondensate, ic->Seed,
		 ic->BumpX, ic->BumpY, ic->BumpRadius, ic->BumpHeight);

  else {

    if (myid == 0) printf("set_initial_conditions: unspecified i.c.\n");
    MPI_Finalize();
    exit(1);
  }

  free(A);

}


/*-----------------------------------------------------------*/

void Tophat(fftw_complex **A, double x0, double y0, double r, double h)
{
  int i1, i2, j1, j2;
  int i, j;
  double x,y;
  double offset = 0;        /* 0.5 for cell-center */

  i1 = (x0 - r)/dx + gp - 1;
  i2 = (x0 + r)/dx + gp + 1;
  j1 = (y0 - r)/dx - myid*ny + gp - 1;
  j2 = (y0 + r)/dx - myid*ny + gp + 1;

  if (i1 < 0)        i1=0;
  if (i2 > nx+2*gp)  i2=nx+2*gp;
  if (j1 < 0)        j1=0;
  if (j2 > ny+2*gp)  j2=ny+2*gp;

  for (j=j1; j<j2; j++) for (i=i1; i<i2; i++) {

    x = (i - gp + offset)*dx;
    y = (myid*ny + j - gp + offset)*dx;

    if ( ((x-x0)*(x-x0) + (y-y0)*(y-y0)) < r*r ) A[j][i][0] = h;
  }
}
/*-----------------------------------------------------------*/

void Gauss(fftw_complex **A, double x0, double y0, double r, double h, double phase)
{
  superGauss(A, x0, y0, r, h, 2, phase);
}

/*-----------------------------------------------------------*/

void superGauss(fftw_complex **A, double x0, double y0, double r, double h, double p, double phase)
{
  int i, j;
  double x,y, rr, rr0;
  double rx, ry, lx, ly;
  double a, cos_ph, sin_ph;
  double offset = 0;        /* 0.5 for cell-center */

  if (r<=0) return;

  cos_ph=cos(phase);
  sin_ph=sin(phase);

  rr0 = r*r;

  for (j=0; j<ny+gp; j++) for (i=0; i<nx+gp; i++) {

    x = (i - gp + offset)*dx;
    y = (myid*ny + j - gp + offset)*dx;

    rx = fabs(x-x0);
    ry = fabs(y-y0);

    if (rx > 0.5*Lx) rx = Lx-rx;  
    if (ry > 0.5*Ly) ry = Ly-ry;  

    rr = rx*rx + ry*ry;
    a = exp ( -pow( rr/rr0, 0.5*p) ); 
    A[j][i][0] += h*a*cos_ph;
    A[j][i][1] += h*a*sin_ph;

  }

}

/*-----------------------------------------------------------*/

void MakeHole(fftw_complex **A, double x0, double y0, double R)
{
  int i, j;
  double x,y, rr, r6, R6, f;
  double rx, ry, lx, ly;
  double offset = 0;        /* 0.5 for cell-center */

  if (R<=0) return;

  R6 = pow(R,6);

  for (j=0; j<ny+gp; j++) for (i=0; i<nx+gp; i++) {

    x = (i - gp + offset)*dx;
    y = (myid*ny + j - gp + offset)*dx;

    rx = fabs(x-x0);
    ry = fabs(y-y0);

    if (rx > 0.5*Lx) rx = Lx-rx;  
    if (ry > 0.5*Ly) ry = Ly-ry;  

    rr = rx*rx + ry*ry;
    r6 = rr*rr*rr;

    f  = 1 - exp(-r6/R6);

    A[j][i][0] = A[j][i][0] * f;
    A[j][i][1] = A[j][i][1] * f;

  }

}


/*-----------------------------------------------------------*/

double HoleRadius(double r, double h, double b)
{
  double R, q, f;

  if (b<=0) return(0.);

  f = (sqrt(14) - sqrt(2))/4;

  q = pow( -log(1-f), -1/6.);

  R = r * q * sqrt( log(2*h/b) );
  
  return(R);

}

/*-----------------------------------------------------------*/
void icGaussNoise(fftw_complex **A, double kmin, double kmax, 
		  double numParticles,  double nCondensate,  int seed,
                  double BumpX, double BumpY, 
		  double BumpRadius, double BumpHeight) {

  double r;

  if (numParticles > 0) {

    icSpectrum(A, kmin, kmax, numParticles, nCondensate, seed, SpGauss);

    r = HoleRadius(BumpRadius, BumpHeight, sqrt(numParticles)/Lx );

    MakeHole(A, BumpX, BumpY, r);

  }

  Gauss(A, BumpX, BumpY, BumpRadius, BumpHeight, 0);

}


/*-----------------------------------------------------------*/

void Ellipse(fftw_complex **A, double x0, double y0, double a, double b, double h)
{
  int i, j;
  double x,y, rr;
  double rx, ry, lx, ly;
  double offset = 0;        /* 0.5 for cell-center */

  for (j=0; j<ny+gp; j++) for (i=0; i<nx+gp; i++) {

    x = (i - gp + offset)*dx;
    y = (myid*ny + j - gp + offset)*dx;

    rx = fabs(x-x0);
    ry = fabs(y-y0);

    if (rx > 0.5*Lx) rx = Lx-rx;  
    if (ry > 0.5*Ly) ry = Ly-ry;  

    rr = rx*rx/(a*a) + ry*ry/(b*b); 
    A[j][i][0] += h*exp(-rr);
  }

}
/*-----------------------------------------------------------*/

void Vortex(fftw_complex **A, double x0, double y0, double r0, double h)
{
  int i, j;
  double x, y, r;
  double cosphi, sinphi, amp;
  double offset = 0;        /* 0.5 for cell-center */

  for (j=0; j<ny+gp; j++) for (i=0; i<nx+gp; i++) {

    x = (i - gp + offset)*dx;
    y = (myid*ny + j - gp + offset)*dx;

    r = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));

    if ((r<=r0) && (r != 0)) {

      cosphi = (x-x0)/r;
      sinphi = (y-y0)/r;
      amp    = h*r/r0*(1-r/r0)*4;

      A[j][i][0] = amp*cosphi;
      A[j][i][1] = amp*sinphi;
    }

  }
}

/*-----------------------------------------------------------*/

void MultiGauss(fftw_complex **A, double rmax, double hmax, int nbumps, int seed)
{
  int k;
  double r, h, x, y;
  double rmin=2;

  srand48(seed);

  for (k=0; k<nbumps; k++){

    x = Lx * drand48();
    y = Ly * drand48();
    r = rmin + (rmax-rmin) * drand48();
    h = hmax * (2*drand48() - 1);

    Gauss(A, x, y, r, h, 0);
  }
}


/*-----------------------------------------------------------*/

void GaussArray(fftw_complex **A, double numParticles, double r, int nbumps, int seed)
{
  double        x, y, pi, h, phase;
  int           i, j, n, N;

  N = nx;
  n = ny;

  srand48(seed);
  pi = acos(0)*2;
  h  = 1.0;

  for (i=0; i<nbumps; i++){
    for (j=0; j<nbumps; j++){

      phase = 2*pi*drand48();

      x = Lx * (i + 0.5)/nbumps;
      y = Ly * (j + 0.5)/nbumps;

      Gauss(A, x, y, r, h, phase);
    }
  }

  normalizeField(A, numParticles);

}

/*-----------------------------------------------------------*/

void GaussHex(fftw_complex **A, double numParticles, double d, double c, int m, int seed)
{
  double        x, y, pi, h, phase;
  double        a, b;
  double        x0, y0, rmax;
  int           i, j, n, N;

  N = nx;
  n = ny;

  srand48(seed);
  pi = acos(0)*2;
  h  = 1.0;

  a=d/2;
  b=d*sqrt(3)/2;

  rmax = (m -1 + 1.e-6)*d;

  for (j=-m; j<m+1; j++){
    for (i=-m; i<=m+1; i++){

      x = i*d + d/2 *(2*(j/2) - j);
      y = j*b; 

      if ( x*x + y*y <= rmax*rmax ) {
	phase = 2*pi*drand48();
	Gauss(A, Lx/2+x, Ly/2+y, c*d, h, phase);
      }
    }
  }

  normalizeField(A, numParticles);

}

/*-----------------------------------------------------------*/

void superGaussHex(fftw_complex **A, double numParticles, double d, double c, double p, int m, 
                    double chirp, double width, double phaserand, double seed)
{
  double        x, y, pi, h, phase;
  double        a, b, rr;
  double        x0, y0, rmax;
  int           i, j, n, N;

  N = nx;
  n = ny;

  srand48(seed);
  pi = acos(0)*2;

  a=d/2;
  b=d*sqrt(3)/2;

  rmax = (m -1 + 1.e-6)*d;

  for (j=-m; j<m+1; j++){
    for (i=-m; i<=m+1; i++){

      x = i*d + d/2 *(2*(j/2) - j);
      y = j*b; 
      rr = x*x + y*y;
      h  = exp( -rr/(width*width));

      if ( rr <= rmax*rmax ) {
	phase = rr * chirp;
        phase = phase + pi*phaserand*(2*drand48() - 1);
	superGauss(A, Lx/2+x, Ly/2+y, c*d, h, p, phase);
      }
    }
  }

  normalizeField(A, numParticles);

}




/*-----------------------------------------------------------*/

void normalizeField(fftw_complex **A, double numParticles)
{
  double        a, b, hh, hhtot, w;
  int           i, j, n, N;
  fftw_complex *data;

  data = A[0];

  N = nx;
  n = ny;

  /*-- normalize the field --*/

  hhtot = 0;

  for (i=0; i<N*n; i++) {
    a = data[i][0];
    b = data[i][1];
    hhtot += a*a + b*b;
  }

  hhtot = hhtot*Lx/N*Ly/N;

  MPI_Allreduce(&hhtot, &hh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  w = sqrt(numParticles/hh);

  for (i=0; i<N*n; i++) {
    data[i][0] = w*data[i][0];
    data[i][1] = w*data[i][1];
  }

}




/*-----------------------------------------------------------*/

void icSpectrum(fftw_complex **A, double kmin, double kmax, 
		double numParticles, double nCondensate, int seed, int SpType)
{
  double        mu, pi, k, kk, hh, h, p, hhtot, w;
  int           N, n, i, j, kx, ky;
  fftw_complex *data;

  data = A[0];

  N = nx;
  n = ny;


  /*-- shape the spectrum --*/

  srand48(seed);
  mu = kmin*kmin;
  pi = 2*asin(1);

  p = 0;

  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    kk = kx*kx + ky*ky;
    k  = sqrt(kk);

    if ( (kmax > 0) && ( k >= kmax ) ) {
      data[N*j+i][0] = 0;
      data[N*j+i][1] = 0;
      continue;
    }

    switch(SpType){
    case SpBump:
      if (k>kmin)   
	h  = (k-kmin)*(k-kmax);
      else 
	h = 0;
      hh = h*h;
      break;
    case SpTail:
      if (k>kmin)    
        hh = 1/kk;   // hh = 1/(kk+mu);
      else
        hh = 0;
      h  = sqrt(hh); 
      break;
    case SpExp:
      h  = exp(-k/kmin);
      hh = h*h;
      break;
    case SpGauss:
      h  = drand48();
      h  = sqrt(-2*log(h));
      h  = h*exp(-kk/mu);
      hh = h*h;
      break;
    }

    p = drand48()*2*pi;

    hhtot += hh;

    data[N*j+i][0] = h*cos(p);
    data[N*j+i][1] = h*sin(p);

  }

  hhtot = hhtot*Lx/N*Ly/N;

  /*-- normalize the spectrum --*/

  MPI_Allreduce(&hhtot, &hh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  w = sqrt(numParticles/hh);

  for (i=0; i<N*n; i++) {

    data[i][0] = w*data[i][0];
    data[i][1] = w*data[i][1];

  }

  /*-- add condensate --*/

  p = drand48()*2*pi;

  if (myid == 0){
    data[0][0] += sqrt(nCondensate)*cos(p)*N;
    data[0][1] += sqrt(nCondensate)*sin(p)*N;
  }

  fft_wrap(data, grid, BACKWARD);

}

/*-----------------------------------------------------------*/
/*-----------------------------------------------------------*/
