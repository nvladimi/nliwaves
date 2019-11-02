
/*---------------------------------------------------------------*/
/*                                                               */
/*                                                               */
/* This code models isolated radially symmetric Keller-Segel     */
/* collapses in the radial coordinates.  It solves the equation  */
/*                                                               */
/*              du/dt = u_rr + (u-1)/r u_r                       */
/*                                                               */
/* with I.C.:     u(0) = 4A r^2/(1+r^2)                          */
/*                                                               */
/* using Runger-Kutta time integration of 4rd order. The spatial */
/* derivatives are computed spectrally.  The boundary conditions */
/* are perioic in [0,L] domain. The code is stable at CFL = 0.28 */
/*                                                               */
/*  Compiling on Pequena:                                        */
/*                                                               */
/*    module load  fftw/3.2.1/intel                              */
/*    icc -O3 kelseg.c -lfftw3 -lm -o kelseg.x                   */
/*                                                               */
/*  Compiling on Mallorn:                                        */
/*                                                               */
/*    gcc -O3 kelseg.c -I/usr/local/fftw/fftw-3.2.2/include      */
/*        -L/usr/local/fftw/fftw-3.2.2/lib -lfftw3 -lm           */
/*                                                               */
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
static void   init_mass();
static void   print_header();

static void   compute_rhs(int N);
static void   euler(int N, double dt);

static double runge_kutta(int N, double t, double dt,
			  int save_center, int save_profile);

static double diagnostics(int N, double t, 
			  int save_center, int save_profile);

static void   rk_zero_out(int N);
static void   rk_add_slope(int N, int w);
static void   rk_part_step(int N, double dt);
static void   rk_full_step(int N, double dt);

static void   reset_fft(int N);
static void   refine_grid(int N);
static void   reset_coord(int N);


/*------------------------------------------------------------------*/

/*-- input parameters and their defaults --*/

double   L        = 25.6;         // length of the box
double   A        =  1.0;         // initial height, A=rho/8
double   tend     = 1000;         // end time
double   cfl      =  0.2;         // specrtal derivatives + RK4
double   ppl      =    4;         // points per L, refinement criterium

int      N0       =  256;         // points at coarsest grid
int      Ngrids   =    3;         // number of grids
int      nout1    =   50;         // frequency of output (center)
int      nout2    = 5000;         // frequency of output (profile) 


/*-- arrays of run-time variables  --*/

char          fbase[80];          // name base for output files

double        *R;                 // coordinate
double        *Uo;                // solution stored
double        *U;                 // RK: intermediate solution
double        *dU;                // RK: avg slope (RHS) collected 


/*--  arrays of FFT variables --*/

fftw_complex  *Psi, *D1, *D2;
fftw_plan      fftPsi, ifftD1, ifftD2, ifftPsi;

int           nio  = 0;            // output file number


/*------------------------------------------------------------------*/

int main(int argc, char **argv)
{

  int          n    = 0;            // time step
  int          grid = 0;            // grid number
  double       t    = 0;            // time   
  double       r    = 1;            // radius, "ell"

  int          N;                   // grid size
  double       dt;                  // time step
  double       dx;                  // resolution

  int          save_center;
  int          save_profile;

  parameters(argc, argv);
  print_header();
  allocate_arrays();
  init_mass();


  N  = N0;
  dx = L/N;
  dt = cfl*dx*dx;

  reset_fft(N);

  /*-- main loop --*/

  while (t <= tend) {


    /*-- update solution and save data after first fft --*/

    if (remainder(n, nout1) == 0) save_center = 1;
    else save_center = 0;

    if (remainder(n, nout2) == 0) save_profile = 1;
    else save_profile = 0;

    r = runge_kutta(N, t, dt, save_center, save_profile);

    n++;
    t+=dt;


    /*-- refine grid if needed --*/

    if ( r/dx < ppl ) {

      if (grid+1 == Ngrids ) {
	printf("Need grid N= %6d at t= %18.12e  ell = %18.12e\n", 2*N, t, r);
	printf("Exiting\n");
	exit(0);
      }

      refine_grid(N);
      N  = N*2;
      dx = dx/2;
      dt = dt/4;
      reset_coord(N);
      grid++;
      printf("regrid to N= %6d at t= %18.12e  ell = %18.12e\n", N, t, r);
 
    }

  }

  exit(0);

}


/*---------------------------------------------------------------*/
void parameters(int argc, char **argv)
{
  int i;

  strcpy (fbase, "ks");

  for( i=0; i<argc; i++){
    if( strcmp(argv[i],"-f") == 0 )           strcpy(fbase,argv[i+1]);

    if( strcmp(argv[i],"-A") == 0 )           A =       atof(argv[i+1]);
    if( strcmp(argv[i],"-L") == 0 )           L =       atof(argv[i+1]);
    if( strcmp(argv[i],"-N0") == 0 )          N0 =      atoi(argv[i+1]);
    if( strcmp(argv[i],"-Ngrids") == 0 )      Ngrids =  atoi(argv[i+1]);

    if( strcmp(argv[i],"-tend") == 0 )        tend =    atof(argv[i+1]);
    if( strcmp(argv[i],"-nout1") == 0 )       nout1 =   atoi(argv[i+1]);
    if( strcmp(argv[i],"-nout2") == 0 )       nout2 =   atoi(argv[i+1]);

    if( strcmp(argv[i],"-cfl") == 0 )         cfl =     atof(argv[i+1]);
    if( strcmp(argv[i],"-ppl") == 0 )         ppl =     atof(argv[i+1]);
  }

    
  printf("\n  started with parameters:\n\n");

  printf("   %s\n", argv[0]);
  printf("      -f %s  -A %f  -L %f  -N0 %d  -Ngrids %d\n", 
	 fbase, A, L, N0, Ngrids);
  printf("      -ppl %4.1f  -cfl %4.2f  -nout1 %d  -nout2 %d  -tend %6.1f\n\n", 
	 ppl, cfl, nout1, nout2, tend);


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
  strcat(fname, ".out");

  fid = fopen (fname, "w");

  fprintf(fid, "%%\n%% Run \"%s\"  started at %s", fbase, date);
  fprintf(fid, "%%\n%% Collapse parameters\n");
  fprintf(fid, "%%      A      = %f\n", A);
  fprintf(fid, "%%      tend   = %f\n", tend);
  fprintf(fid, "%%\n%% Discretization\n");
  fprintf(fid, "%%      L      = %f\n", L);
  fprintf(fid, "%%      N0     = %d\n", N0);
  fprintf(fid, "%%      Ngrids = %d\n", Ngrids);
  fprintf(fid, "%%      cfl    = %f\n", cfl);
  fprintf(fid, "%%      ppl    = %f\n", ppl);
  fprintf(fid, "%%\n%% Frequency of output\n");
  fprintf(fid, "%%      nout1  = %d\n", nout1);
  fprintf(fid, "%%      nout2  = %d\n", nout2);
  fprintf(fid, "%%\n");

  fprintf(fid, "%%_1.t  2.ell  3.dens(0)  4.mass(0)\n");
 
  fclose(fid);
}


/*---------------------------------------------------------------*/

void allocate_arrays()
{

  int Nmax;

  Nmax =  N0 * pow(2, Ngrids-1);

  R    = (double*)malloc( Nmax * sizeof(double));
  Uo   = (double*)malloc( Nmax * sizeof(double));
  U    = (double*)malloc( Nmax * sizeof(double));
  dU   = (double*)malloc( Nmax * sizeof(double));

  D1   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
  D2   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
  Psi  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);

}

/*---------------------------------------------------------------*/

void init_mass()
{

  int i, k;
  double r, rr, dx;

  dx = L/N0;

  for (i=0; i<N0; i++){
    if (i<=N0/2)  k=i; else k=i-N0; 
    r     = k*dx;
    rr    = r*r;
    R[i]  = r;
    Uo[i] = A*4*rr/(1+rr);
  }

}

/*---------------------------------------------------------------*/

void reset_coord(int N)
{

  int i, k;
  double dx;

  dx = L/N;

  for (i=0; i<N; i++){
    if (i<=N/2)  k=i; else k=i-N; 
    R[i]  = k*dx;
  }

}


/*---------------------------------------------------------------*/

void reset_fft(int N)
{

  fftPsi = fftw_plan_dft_1d(N, Psi, D2, FFTW_FORWARD,  FFTW_ESTIMATE);
  ifftD1 = fftw_plan_dft_1d(N, D1,  D1, FFTW_BACKWARD, FFTW_ESTIMATE);
  ifftD2 = fftw_plan_dft_1d(N, D2,  D2, FFTW_BACKWARD, FFTW_ESTIMATE);

}


/*--------------------------------------------------------------*/

void refine_grid(int N)
{

  int  i;

  ifftPsi = fftw_plan_dft_1d(2*N, D2, Psi, FFTW_BACKWARD,  FFTW_ESTIMATE);

  /*-- fill in complex array --*/

  for (i=0; i<N; i++){
    Psi[i][0] = Uo[i];
    Psi[i][1] = 0;
  }


  /*-- double the number of modes --*/

  fftw_execute(fftPsi);

  for (i=N/2+1; i<N; i++) {
    D2[i+N][0] = D2[i][0];
    D2[i+N][1] = D2[i][1];
  }

  for (i=N/2+1; i<N/2+N; i++) {
    D2[i][0] = 0;
    D2[i][1] = 0;
  }

  fftw_execute(ifftPsi);


  /*-- fill in U array --*/

  for (i=0; i<2*N; i++){
    Uo[i] = Psi[i][0]/N;
   }

  reset_fft(2*N);

}

/*--------------------------------------------------------------*/


void rk_zero_out(int N)           // U=0,  dU=0;
{
  int i;
  for (i=0; i<N; i++) { 
     U[i]   = 0; 
    dU[i]   = 0; 
  } 
}


void rk_add_slope(int N, int w)     // dU = dU + U*h,   U holds RHS
{
  int i;
  for (i=0; i<N; i++) 
    dU[i] += U[i]*w; 
}


void rk_part_step(int N, double dt)  // U = Uo + U*h,   RHS is replaced by solution
{
  int i;
  for (i=0; i<N; i++)
    U[i] = Uo[i] + U[i]*dt;
}


void rk_full_step(int N, double dt)  // Uo = Uo + dU*h, 
{
  int i;
  for (i=0; i<N; i++)
    Uo[i] += dU[i]*dt;
}

/*--------------------------------------------------------------*/

void euler(int N, double dt)   // Uo := Uo + dU
{
  int i;

  rk_zero_out(N);

  rk_part_step(N, 0);
  compute_rhs(N);
  rk_add_slope(N, 1);

  rk_full_step(N, dt);


}


/*--------------------------------------------------------------*/

double runge_kutta(int N, double t, double dt, 
		   int save_center, int save_profile)
{

  double ell;

  rk_zero_out(N);

  /*-- stage 1 --*/

  rk_part_step(N, 0);
  compute_rhs(N);

  ell = diagnostics(N, t, save_center, save_profile);

  rk_add_slope(N, 1);


  /*-- stage 2 --*/

  rk_part_step(N, dt/2);
  compute_rhs(N);
  rk_add_slope(N, 2);


  /*-- stage 3 --*/

  rk_part_step(N, dt/2);
  compute_rhs(N);
  rk_add_slope(N, 2);


  /*-- stage 4 --*/

  rk_part_step(N, dt);
  compute_rhs(N);
  rk_add_slope(N, 1);


 /*-- update solution --*/

  rk_full_step(N, dt/6);
  Uo[0] = 0;

  return(ell);

}


/*---------------------------------------------------------------*/

void compute_rhs(int N)   // rewrites U by rhs(U)
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
    Psi[i][1] = 0;
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


  /*-- assemble RHS --*/

  for (i=1; i<N; i++)
    U[i]   = D2[i][0] + (U[i] - 1) * D1[i][0] / R[i];
  
}

/*--------------------------------------------------------------*/

double diagnostics(int N, double t, int save_center, int save_profile){

  FILE *fid;
  char   fname[80];
  int    i, k;

  double  rho, ell, mass;

  /*-- find density at the center and mass at the border --*/

  rho = D2[0][0];

  ell = sqrt(8/rho);

  mass = Uo[N/2+1];


  /*-- write density --*/ 

  if (save_center) {

    strcpy(fname, fbase);
    strcat(fname, ".out");

    fid = fopen (fname, "at");
    fprintf(fid,  "  %18.12e  %18.12e  %18.12e %18.12e\n",  
	    t, ell, rho, mass);
    fclose(fid);

  }


  /*-- write profile information --*/

  if (save_profile) {

    sprintf(fname,"%s.%04d", fbase, nio);
    
    fid = fopen (fname, "w");

    fprintf(fid, "%%1.r  2.rho  3.m  4.m_r  5.m_rr\n" );
    fprintf(fid, "%%time = %12.8f  ell = %12.8f  mass = %12.8f\n\n", 
	    t, ell, mass);

    fprintf(fid, "%16.8e  %16.8e  %16.8e  %16.8e  %16.8e\n", 
	    R[0], rho, Uo[0], D1[0][0], D2[0][0]);

    for (i=1; i<=N/2; i++)  
      fprintf(fid, "%16.8e  %16.8e  %16.8e  %16.8e  %16.8e\n", 
	      R[i], D1[i][0]/R[i], Uo[i], D1[i][0], D2[i][0]);

    fclose(fid);

    nio++;
  }


  return(ell);
 
}


/*--------------------------------------------------------------*/
