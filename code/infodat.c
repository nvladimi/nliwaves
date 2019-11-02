#include "header.h"

static  int 		 myid, np;
static  int 		 N0;
static  double           L;
static  double           pi;

static  fftw_complex    *data;
static  fftw_complex    *fdata;

static  double           vortexTh;
static  double           clpsPsi0;
static  int              collapses;        
static  char             filename[80];
static  char             msg[80];


double astigmatism(int grid);

/*----------------------------------------------------------------*/

void info_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag,
	       fftw_complex *psi, fftw_complex *psihat){

  FILE      *thefile;
  int k;

  data  = psi;
  fdata = psihat;

  /*-- set parameters --*/

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  N0 =  geom->N0;
  L  =  geom->Lx;
  pi =  acos(0)*2;

  vortexTh = diag->vortexTh;
  collapses = diag->clpsCmax;
  clpsPsi0  = diag->clpsPsi0;

  strcpy(filename, ctrl->runname);
  strcat(filename, ".dat");

  /*-- print header --*/

  if (myid == 0) {
    thefile = fopen (filename, "a");
    fprintf(thefile,  "%%\n%% 1.time  2.|psi|_max  3.|psi|_avg  4.N  5.E_kin");
    fprintf(thefile,  "  6.(-focus)*E_pot  7.condensate  8.astigmatism \n%%\n");
    fclose(thefile); 
  }
}


/*-----------------------------------------------------------*/

void info_output(int grid, double t)
{
  FILE          *thefile;

  double         psi, psiavg, psisq, psiquad;
  double         No, vortices, Ek, Ep, A;
  double         u, v, sq, max;

  double         a        = 0;
  double         b        = 0;
  double         mx       = 0;
  double         sum      = 0;
  double         sumsq    = 0;
  double         sumquad  = 0;

  int            i, n, nv, N;

  N = N0*pow(2,grid);
  n = N*N/np;

  /*-- local sums --*/ 

  for (i=0; i<n; i++) {

      u = data[i][0];
      v = data[i][1];
      sq = u*u + v*v;

      a       += u;
      b       += v;
      sum     += sqrt(sq);
      sumsq   += sq;
      sumquad += sq*sq;

      if (sq>mx) mx=sq;

  }

  mx = sqrt(mx);

  /*-- global sums --*/

  MPI_Reduce(&a,       &u,       1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&b,       &v,       1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum,     &psi,     1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumsq,   &psisq,   1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumquad, &psiquad, 1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Allreduce(&mx,   &max,     1,  MPI_DOUBLE, MPI_MAX,    MPI_COMM_WORLD);


  /*-- rescale to represent total or average quantities --*/

  psiavg  =  psi / (N*N);
  psisq   =  psisq * (L*L)/(N*N);
  psiquad =  psiquad * (L*L)/(N*N);
  u       =  u / (N*N);
  v       =  v / (N*N);
  No      =  u*u+v*v;

  /*-- vortex fraction --*/

  if (vortexTh) {

    mx = vortexTh*psiavg;
    mx = mx*mx;
    nv  = 0;

    for (i=0; i<n; i++) {

      u = data[i][0];
      u = data[i][1];
      if ( u*u + v*v  < mx ) nv++ ;

    }  

    a = ((double)nv)/(N*N);
    MPI_Reduce(&a, &vortices, 1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  } else vortices = 0;


  /*-- output to file --*/

  Ek = Ekin(grid);
  Ep = 0.5*psiquad;

  if (collapses == 1) A = astigmatism(grid); else A=0;

  if (myid == 0) {

    thefile = fopen (filename, "a");
    fprintf(thefile, "%18.12e  %18.12e  %18.12e  %18.12e  ", 
	    t, max, psiavg, psisq);  
    fprintf(thefile, "%18.12e  %18.12e  %18.12e  %18.12e\n", 
	    Ek, Ep, No, A);  
    fclose(thefile); 

  }


#ifdef EXIT_AT_PSIMAX
  if (max > clpsPsi0) {
      sprintf(msg, "#msg: psi maximum reached at t = %18.12e\n", t);
      io_message(msg,0);  MPI_Finalize();  exit(0);
  }
#endif


}

/*----------------------------------------------------------------*/

double abspsi_max(int grid)
{
  int            i, n, N;
  double         a, b, sq;
  double         max=0;

  N = N0*pow(2,grid);
  n = N*N/np;

  for (i=0; i<n; i++) {

      a = data[i][0];
      b = data[i][1];
      sq = a*a + b*b;

      if (sq>max) max=sq;
  }

  MPI_Reduce(&max, &sq, 1,  MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  return(sqrt(sq));

}

/*----------------------------------------------------------------*/

double Ekin(int grid)
{
  double        c, u, v, s;
  int           N, n, i, j, kx, ky;

  N = N0*pow(2,grid);
  n = N/np;

  c = 4*pi*pi/(N*N);
  s = 0;

  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    u = fdata[N*j+i][0];
    v = fdata[N*j+i][1];

    s += (u*u + v*v)*(kx*kx + ky*ky);

  }

  s = s*c;

  MPI_Reduce(&s, &c,  1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return(c);

}

/*-----------------------------------------------------------*/


double astigmatism(int grid)
{
  double        c, u, v, sum_x, sum_y;
  int           N, n, i;

  N = N0*pow(2,grid);
  n = N/np;

  sum_x = 0;
  sum_y = 0;

  /*-- x-derivative --*/

  for (i=0; i<N*n; i++) {
    u = data[i][0];
    v = data[i][1];
    fdata[i][0] = u*u + v*v;
    fdata[i][1] = 0;
  }

  fft_deriv_x(fdata, grid);

  for (i=0; i<N*n; i++) {
    u = fdata[i][0];
    v = fdata[i][1];
    sum_x += sqrt(u*u + v*v);
  }

  /*-- y-derivative --*/

  for (i=0; i<N*n; i++) {
    u = data[i][0];
    v = data[i][1];
    fdata[i][0] = u*u + v*v;
    fdata[i][1] = 0;
  }

  fft_deriv_y(fdata, grid);

  for (i=0; i<N*n; i++) {
    u = fdata[i][0];
    v = fdata[i][1];
    sum_y += sqrt(u*u + v*v);
  }


  /*-- ratio --*/

  MPI_Reduce(&sum_x, &c,  1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  sum_x = c;

  MPI_Reduce(&sum_y, &c,  1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  sum_y = c;
 
  c = sum_x / sum_y;

  fft_make_fourier(data, fdata, grid);

  return(c);

}

/*----------------------------------------------------------------*/
