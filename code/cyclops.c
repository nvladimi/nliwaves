
/*---------------------------------------------------------------*/
/*                                                               */
/*                                                               */
/* This code models isolated radially symmetric collapses in     */
/* the radial coordinates.  It solves the focusing NLS equation  */
/*                                                               */
/*      i dPsi/dt + (1-i*a*eps)*del^2(Psi) +                     */
/*                                                               */
/*           + (1+i*c*eps*|Psi|^s)*|Psi|^2*Psi = i*b*eps*Psi     */
/*                                                               */
/* using Runger-Kutta time integration of 4rd order. The spatial */
/* derivatives are computed spectrally.  The boundary conditions */
/* are perioic in [-L/2, L/2] domain. The initial conditions are */
/* Gaussian for Re(Psi) and zero for Im(Psi).  The width, r0,    */
/* and the height, h0, of the initial Gaussian are input         */
/* parameters, as well the coefficients in the a,b,c,s, and eps  */
/* in the equation.                                              */
/*                                                               */
/* The code is stable at CFL = 0.28, and dx=1/160 is enough to   */
/* resolve collapses up to |psi|_max = 50.                       */
/*                                                               */
/*  Compiling on Pequena:                                        */
/*                                                               */
/*    module load  fftw/3.2.1/intel                              */
/*    icc -O3 cyclops.c -lfftw3 -lm -o cyclops.x                 */
/*                                                               */
/*  Compiling on Mallorn:                                        */
/*                                                               */
/*    gcc -O3 cyclops.c -I/usr/local/fftw/fftw-3.2.2/include     */
/*        -L/usr/local/fftw/fftw-3.2.2/lib -lfftw3 -lm           */
/*                                                               */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>

static void   parameters(int argc, char **argv);

static void   allocate_arrays();
static void   init_psi();
static void   read_psi();
static void   print_header();

static void   compute_laplacian(int N);            // ddU = laplace(U)
static void   compute_rhs(int N);                  // U   = RHS(U, ddU)
static void   runge_kutta(int N, double dt);       // Uo  = Uo + dU
static void   euler(int N, double dt);             // Uo  = Uo + dU

static void   rk_zero_out(int N);
static void   rk_add_slope(int N, int w);
static void   rk_part_step(int N, double dt);
static void   rk_full_step(int N, double dt);

static double output_center(double t);
static void   output_profile(int N, int nio, double t);
static void   beam_quality(int N, double t, double h);


static void   reset_fft(int N);
static void   refine_grid(int N);
static void   derefine_grid(int N);
static void   reset_coord(int N);
static void   save_center(int N);

/*------------------------------------------------------------------*/

/*-- input parameters and their defaults --*/

double   a        =  1;           // coefficient "a" in equation 
double   b        =  0;           // coefficient "b" in equation 
double   c        =  1;           // coefficient "c" in equation 
double   s        =  0;           // exponent    "s" in equation 
double   eps      =  0.01;        // coefficient "epsilon"

double   L        = 12.8;         // length of the box
double   r0       =  1.0;         // initial width of gaussian psi
double   h0       =  3.2;         // initial height of gaussian psi
double   tend     = 10.0;         // end time
double   toff     = 80.0;         // time to turn off nonlinearity
double   ton      = 90.0;         // time to turn on  nonlinearity
double   cfl      =  0.2;         // specrtal derivatives + RK4
double   hdx      =  0.3;         // refinement criterium

int      N0       =  256;         // points at coarsest grid
int      Ngrids   =    3;         // number of grids
int      nout1    =   32;         // frequency of output (center)
int      nout2    =  320;         // frequency of output (profile) 
int      linear   =  0;           // switch to linear propagation
int      use_ic   =  0;           // use IC from file


/*-- arrays of run-time variables  --*/

char          fbase[80];          // name base for output files

double        *R;                 // coordinate
double        *Uo, *Vo;           // solution stored

double        *ddU, *ddV;         // RK: laplacian in RHS
double        *dU,  *dV;          // RK: avg slope (RHS) collected 
double        *U,   *V;           // RK: intermediate solution

double        ddUc[6], ddVc[6];   // psi_rr at 6 center points
double          Uc[6],   Vc[6];   // psi    at 6 center points


/*--  arrays of FFT variables --*/

fftw_complex  *Psi, *D1, *D2;
fftw_plan      fftPsi, ifftD1, ifftD2, ifftPsi;


/*------------------------------------------------------------------*/

int main(int argc, char **argv)
{

  int          n    = 0;            // time step
  int          nio  = 0;            // output file number
  int          grid = 0;            // grid number
  double       t    = 0;            // time   
  double       tmax = 0;            // time at max                  
  double       hmax = 0;            // height at max                  

  int          N;                   // grid size
  double       dt;                  // time step
  double       dx;                  // resolution
  double       h;                   // height


  parameters(argc, argv);
  print_header();
  allocate_arrays();
  if (use_ic) read_psi(); 
  else init_psi();


  N = N0;
  h = h0;

  dx = L/N;
  dt = cfl*dx*dx;

  reset_fft(N);

  /*-- main loop --*/

  while (t <= tend) {

      if (t>=toff)  linear = 1; else linear = 0;
    //if ((t>=toff) && (t<ton)) linear = 1; else linear = 0;
    //if (fmod(t, toff+ton) - ton > 0) linear = 1; else linear = 0;
    //if (t<0.6) linear=0;

    /*-- save collapse profile before we mangled it --*/

    if (remainder(n, nout2) == 0) output_profile(N, nio++, t);

    /*-- finish, if height drops 10% below maximum --*/

    if (h > hmax) {hmax=h; tmax=t;}
    //if (h < 0.9*hmax) exit(0);

    /*-- adjust grid if needed --*/

    if ( h*dx > hdx ){
      if ( grid+1 == Ngrids ) exit(0);
      refine_grid(N);
      N  = N*2;
      dx = dx/2;
      dt = dt/4;
      reset_coord(N);
      grid++;
    }

    if (( h*dx < 0.45*hdx ) && (grid != 0)) {
      derefine_grid(N);
      N  = N/2;
      dx = dx*2;
      dt = dt*4;
      reset_coord(N);
      grid--;
    }

    /*-- update solution and save center data at first fft --*/

    runge_kutta(N, dt);

    /*-- diagnostics at the center uses earlier spectral data --*/

    if (remainder(n, nout1) == 0) {
      h = output_center(t);
      if (t >= toff) beam_quality(N, t, h);
    }

    n++;
    t+=dt;

  }

  exit(0);

}


/*---------------------------------------------------------------*/
void parameters(int argc, char **argv)
{
  int i;

  strcpy (fbase, "cyl");

  for( i=0; i<argc; i++){
    if( strcmp(argv[i],"-f") == 0 )           strcpy(fbase,argv[i+1]);

    if( strcmp(argv[i],"-a") == 0 )           a =       atof(argv[i+1]);
    if( strcmp(argv[i],"-b") == 0 )           b =       atof(argv[i+1]);
    if( strcmp(argv[i],"-c") == 0 )           c =       atof(argv[i+1]);
    if( strcmp(argv[i],"-s") == 0 )           s =       atof(argv[i+1]);
    if( strcmp(argv[i],"-eps") == 0 )         eps =     atof(argv[i+1]);

    if( strcmp(argv[i],"-r0") == 0 )          r0 =      atof(argv[i+1]);
    if( strcmp(argv[i],"-h0") == 0 )          h0 =      atof(argv[i+1]);
    if( strcmp(argv[i],"-L") == 0 )           L =       atof(argv[i+1]);
    if( strcmp(argv[i],"-N0") == 0 )          N0 =      atoi(argv[i+1]);
    if( strcmp(argv[i],"-Ngrids") == 0 )      Ngrids =  atoi(argv[i+1]);

    if( strcmp(argv[i],"-tend") == 0 )        tend =    atof(argv[i+1]);
    if( strcmp(argv[i],"-toff") == 0 )        toff =    atof(argv[i+1]);
    if( strcmp(argv[i],"-ton") == 0 )         ton  =    atof(argv[i+1]);
    if( strcmp(argv[i],"-nout1") == 0 )       nout1 =   atoi(argv[i+1]);
    if( strcmp(argv[i],"-nout2") == 0 )       nout2 =   atoi(argv[i+1]);

    if( strcmp(argv[i],"-cfl") == 0 )         cfl =     atof(argv[i+1]);
    if( strcmp(argv[i],"-hdx") == 0 )         hdx =     atof(argv[i+1]);
    if( strcmp(argv[i],"-use_ic") == 0 )      use_ic =  atoi(argv[i+1]);
  }

   
  printf("\n  started with parameters:\n\n");

  printf("   %s -f %s -use_ic %d\n", argv[0], fbase, use_ic);
  printf("      -a %f -b %f -c %f -s %f -eps %f\n", a, b, c, s, eps);
  printf("      -r0 %f -h0 %f -tend %f\n", r0, h0, tend);
  printf("      -toff %f  -ton %f\n", toff, ton);
  printf("      -L %f -N0 %d -Ngrids %d -hdx %f\n", L, N0, Ngrids, hdx);
  printf("      -cfl %f -nout1 %d -nout2 %d\n\n", cfl, nout1, nout2);


}

/*---------------------------------------------------------------*/

void print_header()
{
  FILE *fid;
  char fname[80];
  time_t     tmr;
  char*      date;


  /*-- get current time --*/

  tmr = time(NULL);
  date = ctime(&tmr); 

  strcpy(fname, fbase);
  strcat(fname, ".clps");

  fid = fopen (fname, "w");

  fprintf(fid, "%%\n%% Run \"%s\"  started at %s", fbase, date);
  fprintf(fid, "%%\n%% Coefficients in the equation\n");
  fprintf(fid, "%%      a      = %f\n", a);
  fprintf(fid, "%%      b      = %f\n", b);
  fprintf(fid, "%%      c      = %f\n", c);
  fprintf(fid, "%%      s      = %f\n", s);
  fprintf(fid, "%%      eps    = %f\n", eps);
  fprintf(fid, "%%\n%% Collapse parameters\n");
  fprintf(fid, "%%      use_ic = %d\n", use_ic);
  fprintf(fid, "%%      r0     = %f\n", r0);
  fprintf(fid, "%%      h0     = %f\n", h0);
  fprintf(fid, "%%      tend   = %f\n", tend);
  fprintf(fid, "%%      toff   = %f\n", toff);
  fprintf(fid, "%%      ton    = %f\n", ton);
  fprintf(fid, "%%\n%% Discretization\n");
  fprintf(fid, "%%      L      = %f\n", L);
  fprintf(fid, "%%      N0     = %d\n", N0);
  fprintf(fid, "%%      Ngrids = %d\n", Ngrids);
  fprintf(fid, "%%      cfl    = %f\n", cfl);
  fprintf(fid, "%%      hdx    = %f\n", hdx);
  fprintf(fid, "%%\n%% Frequency of output\n");
  fprintf(fid, "%%      nout1  = %d\n", nout1);
  fprintf(fid, "%%      nout2  = %d\n", nout2);
  fprintf(fid, "%%\n");

  fprintf(fid, "%%_1.t   2.|psi|  3.|psi|_rr  4.phase_rr  5.phase\n");
  fprintf(fid, "  NaN NaN NaN NaN NaN  %% for compartibility \n\n");
 
  fclose(fid);
}


/*---------------------------------------------------------------*/

void allocate_arrays()
{

  int Nmax;

  Nmax =  N0 * pow(2, Ngrids-1);

  R   = (double*)malloc( Nmax * sizeof(double));
  Uo  = (double*)malloc( Nmax * sizeof(double));
  Vo  = (double*)malloc( Nmax * sizeof(double));

  U   = (double*)malloc( Nmax * sizeof(double));
  V   = (double*)malloc( Nmax * sizeof(double));
  ddU = (double*)malloc( Nmax * sizeof(double));
  ddV = (double*)malloc( Nmax * sizeof(double));
  dU  = (double*)malloc( Nmax * sizeof(double));
  dV  = (double*)malloc( Nmax * sizeof(double));

  D1  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
  D2  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
  Psi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);

}

/*---------------------------------------------------------------*/

void read_psi()
{

  int     i;
  double  r, dx;
  char    fname[80];
  FILE   *fid;
  
  dx = L/N0;

  for (i=0; i<N0; i++)  R[i] = (i+0.5)*dx - L/2;

  sprintf(fname, "%s.ic", fbase);

  if ( (fid = fopen(fname, "rb")) == NULL ) {
    printf ("\n  File \"%s\" does not exist.\n\n", fname);
    exit(1);
  }else{
    if (fread(Uo, N0*sizeof(double), 1, fid) != 1){
      printf ("\n File \"%s\" has wrong size. \n\n", fname ); 
      exit(1);
      }
    if (fread(Vo, N0*sizeof(double), 1, fid) != 1){
      printf ("\n File \"%s\" has wrong size. \n\n", fname ); 
      exit(1);
      }
    fclose(fid); 
  }

}


/*---------------------------------------------------------------*/

void init_psi()
{

  int i;
  double r, dx;

  dx = L/N0;

  for (i=0; i<N0; i++){
    r = (i+0.5)*dx - L/2;
    R[i]  = r;
    Uo[i] = h0 * exp( -r*r/(r0*r0));
    //Uo[i] = h0 * exp( - pow(r/r0, 8));
    Vo[i] = 0;
  }

}



/*---------------------------------------------------------------*/

void reset_coord(int N)
{

  int i;
  double r, dx;

  dx = L/N;

  for (i=0; i<N; i++){
    r = (i+0.5)*dx - L/2;
    R[i]  = r;
  }

}


/*---------------------------------------------------------------*/

void output_profile(int N, int nio, double t){

  FILE *fid;
  char   fname[80];
  double amp;
  int    i;

  sprintf(fname,"%s.rad.%04d", fbase, nio);

  fid = fopen (fname, "w");

  fprintf(fid, "%%1.r  2.abs(Psi)  3.Re(Psi)  4.Im(Psi)\n" );
  fprintf(fid, "%%time = %12.8f\n\n", t);

  for (i=0; i<N; i++){
    amp = sqrt (Uo[i]*Uo[i] + Vo[i]*Vo[i]);
    fprintf(fid, "%16.8e  %16.8e  %16.8e  %16.8e\n",
	    R[i], amp, Uo[i], Vo[i]);
  }

  fclose(fid);
}

/*---------------------------------------------------------------*/

double output_center(double t){

  FILE *fid;
  char   fname[80];
  double u, v, h, p;
  double ddu, ddv, ddh, ddp;
  int    i0, i1, i2, i3, i4, i5;

  strcpy(fname, fbase);
  strcat(fname, ".clps");


  /*-- Find derivatives --*/

  u   = ( 3*Uc[0] - 25*Uc[1] + 150*Uc[2] + 
          3*Uc[5] - 25*Uc[4] + 150*Uc[3]    ) / 256;
  v   = ( 3*Vc[0] - 25*Vc[1] + 150*Vc[2] + 
          3*Vc[5] - 25*Vc[4] + 150*Vc[3]    ) / 256;

  ddu = ( 3*ddUc[0] - 25*ddUc[1] + 150*ddUc[2] + 
          3*ddUc[5] - 25*ddUc[4] + 150*ddUc[3]    ) / 256;
  ddv = ( 3*ddVc[0] - 25*ddVc[1] + 150*ddVc[2] + 
          3*ddVc[5] - 25*ddVc[4] + 150*ddVc[3]    ) / 256;


  h   = sqrt(u*u + v*v);
  p   = atan2(v,u);
  ddh = (u*ddu + v*ddv)/h;
  ddp = (u*ddv - v*ddu)/(h*h);


  /*-- write data to file --*/

  fid = fopen (fname, "at");
  fprintf(fid,  "  %18.12e %19.12e %19.12e %19.12e %19.12e\n", 
	  t, h, ddh, ddp, p);
  fclose(fid);

  return(h);
 
}


/*---------------------------------------------------------------*/

void beam_quality(int N, double t, double h){

  FILE *fid;
  char   fname[80];

  double  M1,  M2,  M3,  M4;  //  moments of intensity 
  double  N1,  N2,  N3,  N4;  //  optical power

  static double pp1 = 0;      // noise thresholds
  static double pp2 = 0; 
  static double pp3 = 0; 
  static double pp4 = 0; 

  double pi, pp, dn, dm, dx;
  double r, f, pp0;
  int    i;

  strcpy(fname, fbase);
  strcat(fname, ".gbq");

  pp1 = 2.e-2*h*h;
  pp2 = 1.e-2*h*h;
  pp3 = 5.e-3*h*h;
  pp4 = 1.e-3*h*h;

  //if (pp3 == 0)  pp3 =  2.e-2*h*h;
  //if (pp4 == 0)  pp4 =  1.e-2*h*h;

  /*-- beam quality --*/
 
  M1 = 0;  M2 = 0;  M3 = 0;  M4 = 0; 
  N1 = 0;  N2 = 0;  N3 = 0;  N4 = 0; 

  pi = acos(0)*2;
  dx = L/N;

  pp=0;

  for (i=N/2; i<N; i++) {

    pp0 = pp;

    pp  = Uo[i]*Uo[i] + Vo[i]*Vo[i];
    r   = R[i];
 
    dn = pp*r;
    dm = dn*r*r;

    if (pp > pp1) {
      M1 += dm; 
      N1 += dn; 
    } else if (pp0 > pp1) {
      f = (pp1 - 0.5*(pp+pp0))/(pp - pp0);
      M1 += dm*f;
      N1 += dn*f;
     }

    if (pp > pp2) {
      M2 += dm; 
      N2 += dn; 
    } else if (pp0 > pp2) {
      f = (pp2 - 0.5*(pp+pp0))/(pp - pp0);
      M2 += dm*f;
      N2 += dn*f;
    }

    if (pp > pp3) {
      M3 += dm; 
      N3 += dn; 
    } else if (pp0 > pp3) {
      f = (pp3 - 0.5*(pp+pp0))/(pp - pp0);
      M3 += dm*f;
      N3 += dn*f;
     }

    if (pp > pp4) {
      M4 += dm; 
      N4 += dn; 
    } else if (pp0 > pp4) {
      f = (pp4 - 0.5*(pp+pp0))/(pp - pp0);
      M4 += dm*f;
      N4 += dn*f;
    }
    
  }

  /*-- write data to file --*/

  fid = fopen (fname, "at");
  fprintf(fid,  "  %18.12e %19.12e  %12.6e %12.6e  %12.6e %12.6e  %12.6e %12.6e  %12.6e %12.6e \n", 
	  t, h,   
          sqrt(2*M1/N1),   N1*2*pi*dx,
          sqrt(2*M2/N2),   N2*2*pi*dx,
          sqrt(2*M3/N3),   N3*2*pi*dx,
          sqrt(2*M4/N4),   N4*2*pi*dx
         );
  fclose(fid);
 
}

/*--------------------------------------------------------------*/

void reset_fft(int N)
{

  fftPsi = fftw_plan_dft_1d(N, Psi, D2, FFTW_FORWARD,  FFTW_ESTIMATE);
  ifftD1 = fftw_plan_dft_1d(N, D1,  D1, FFTW_BACKWARD, FFTW_ESTIMATE);
  ifftD2 = fftw_plan_dft_1d(N, D2,  D2, FFTW_BACKWARD, FFTW_ESTIMATE);

}


/*--------------------------------------------------------------*/

void refine_grid(int N)
{

  int  i, k;
  double u, v, c, s, pi, q;

  pi = acos(0)*2;
  q = -pi/(2*N);

  ifftPsi = fftw_plan_dft_1d(2*N, D2, Psi, FFTW_BACKWARD,  FFTW_ESTIMATE);

  /*-- fill in complex array --*/

  for (i=0; i<N; i++){
    Psi[i][0] = Uo[i];
    Psi[i][1] = Vo[i];
  }


  /*-- double the number of modes, shift to new cell centers --*/

  fftw_execute(fftPsi);

  for (i=0; i<=N/2; i++) {

    k = i;
    u = D2[i][0];
    v = D2[i][1];
    c = cos(q*k);
    s = sin(q*k);
    D2[i][0] = u*c - v*s;
    D2[i][1] = u*s + v*c;

  }

  for (i=N/2+1; i<N; i++) {

    k = i-N;
    u = D2[i][0];
    v = D2[i][1];
    c = cos(q*k);
    s = sin(q*k);
    D2[i+N][0] = u*c - v*s;
    D2[i+N][1] = u*s + v*c;

  }

  for (i=N/2+1; i<N/2+N; i++) {
    D2[i][0] = 0;
    D2[i][1] = 0;
  }

  fftw_execute(ifftPsi);


  /*-- fill in U and V array --*/

  for (i=0; i<2*N; i++){
    Uo[i] = Psi[i][0]/N;
    Vo[i] = Psi[i][1]/N;
  }

  reset_fft(2*N);

}
/*--------------------------------------------------------------*/

void derefine_grid(int N)
{

  int  i, k;
  double u, v, c, s, pi, q;

  pi = acos(0)*2;
  q = pi/N;

  ifftPsi = fftw_plan_dft_1d(N/2, D2, Psi, FFTW_BACKWARD,  FFTW_ESTIMATE);

  /*-- fill in complex array --*/

  for (i=0; i<N; i++){
    Psi[i][0] = Uo[i];
    Psi[i][1] = Vo[i];
  }


  /*-- double the number of modes, shift to new cell centers --*/

  fftw_execute(fftPsi);

  for (i=0; i<=N/4; i++) {

    k = i;
    u = D2[i][0];
    v = D2[i][1];
    c = cos(q*k);
    s = sin(q*k);
    D2[i][0] = u*c - v*s;
    D2[i][1] = u*s + v*c;

  }

  for (i=3*N/4+1; i<N; i++) {

    k = i-N;
    u = D2[i][0];
    v = D2[i][1];
    c = cos(q*k);
    s = sin(q*k);
    D2[i-N/2][0] = u*c - v*s;
    D2[i-N/2][1] = u*s + v*c;

  }

  fftw_execute(ifftPsi);


  /*-- fill in U and V array --*/

  for (i=0; i<N/2; i++){
    Uo[i] = Psi[i][0]/N;
    Vo[i] = Psi[i][1]/N;
  }

  reset_fft(N/2);

}


/*--------------------------------------------------------------*/

void compute_laplacian(int N)   // ddU = compute_laplacian(U)
{

  double        c0, c1, c2, kc1, kc2, pi;
  int           i, k;

  pi = acos(0)*2;
  c0 = 2*pi/L;
  c1 = c0/N;
  c2 = c0*c0/N;


  /*-- fill in complex array --*/

  for (i=0; i<N; i++){
    Psi[i][0] = U[i];
    Psi[i][1] = V[i];
  }


  /*-- compute derivatives --*/

  fftw_execute(fftPsi);

  for (i=0; i<N; i++) {

    if (i<= N/2) k=i;  else k=i-N;

    kc1 =    k*c1;
    kc2 = -k*k*c2;

    D1[i][0] = - D2[i][1] * kc1;
    D1[i][1] =   D2[i][0] * kc1;

    D2[i][0] =   D2[i][0] * kc2;
    D2[i][1] =   D2[i][1] * kc2;
  }

  fftw_execute(ifftD1);
  fftw_execute(ifftD2);


  /*-- assemble laplacian --*/

  for (i=0; i<N; i++){
    ddU[i] = D2[i][0] + D1[i][0] / R[i];
    ddV[i] = D2[i][1] + D1[i][1] / R[i];
  }

}

/*--------------------------------------------------------------*/

void save_center(int N)  // save psi and psi_rr at center points,
                         // to be called after compute_laplacian 
{
  int i;

  for (i=0; i<6; i++) {
    ddUc[i] = D2[i + N/2 - 3][0];
    ddVc[i] = D2[i + N/2 - 3][1];
    Uc[i]   =  U[i + N/2 - 3];
    Vc[i]   =  V[i + N/2 - 3];
  }

}

/*--------------------------------------------------------------*/

void rk_zero_out(int N)           // U=0,  dU=0;
{
  int i;
  for (i=0; i<N; i++) { 
     U[i]   = 0; 
     V[i]   = 0; 
    dU[i]   = 0; 
    dV[i]   = 0;
  } 
}


void rk_add_slope(int N, int w)     // dU = dU + U*h,   U holds RHS
{
  int i;
  for (i=0; i<N; i++) { 
    dU[i] += U[i]*w; 
    dV[i] += V[i]*w; 
  }
}


void rk_part_step(int N, double dt)  // U = Uo + U*h,   RHS is replaced by solution
{
  int i;
  for (i=0; i<N; i++) { 
    U[i] = Uo[i] + U[i]*dt; 
    V[i] = Vo[i] + V[i]*dt; 
  }
}


void rk_full_step(int N, double dt)  // Uo = Uo + dU*h, 
{
  int i;
  for (i=0; i<N; i++) { 
    Uo[i] += dU[i]*dt; 
    Vo[i] += dV[i]*dt; 
  }
}

/*--------------------------------------------------------------*/

void euler(int N, double dt)   // Uo := Uo + dU
{
  int i;

  rk_zero_out(N);

  rk_part_step(N, 0);
  compute_laplacian(N);
  save_center(N);
  compute_rhs(N);
  rk_add_slope(N, 1);

  rk_full_step(N, dt);


}


/*--------------------------------------------------------------*/

void runge_kutta(int N, double dt)   // Uo := Uo + dU
{

  rk_zero_out(N);

  /*-- stage 1 --*/

  rk_part_step(N, 0);
  compute_laplacian(N);
  save_center(N);
  compute_rhs(N);
  rk_add_slope(N, 1);


  /*-- stage 2 --*/

  rk_part_step(N, dt/2);
  compute_laplacian(N);
  compute_rhs(N);
  rk_add_slope(N, 2);


  /*-- stage 3 --*/

  rk_part_step(N, dt/2);
  compute_laplacian(N);
  compute_rhs(N);
  rk_add_slope(N, 2);


  /*-- stage 4 --*/

  rk_part_step(N, dt);
  compute_laplacian(N);
  compute_rhs(N);
  rk_add_slope(N, 1);


 /*-- update solution --*/

  rk_full_step(N, dt/6);

}

/*--------------------------------------------------------------*/

void compute_rhs(int N)  // U = rhs(U,ddU) 
{

  int i;
  double u, v, ae, be, ce, pp, cePs;


  if (linear) {

    for (i=0; i<N; i++) { 
      U[i] = - ddV[i]; 
      V[i] =   ddU[i];
    }

    return;

  }


  ae = a*eps;
  be = b*eps;
  ce = c*eps;


  for (i=0; i<N; i++){

    u    = U[i];
    v    = V[i];
    pp   = u*u + v*v;

    if (pp>0) cePs = ce*pow(pp, 0.5*s); else cePs=0;  // sePs=|psi|^s

    U[i] = - ddV[i] + ae*ddU[i] + pp*(-v - cePs*u) + be*u;
    V[i] =   ddU[i] + ae*ddV[i] + pp*( u - cePs*v) + be*v;

  }

}
/*---------------------------------------------------------------*/
