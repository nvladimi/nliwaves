#include "header.h"

static  int 		 myid, np;
static  int 		 N0;
static  double           L;
static  double           pi;

static  fftw_complex    *data;
static  fftw_complex    *fdata;

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

 
  strcpy(filename, ctrl->runname);
  strcat(filename, ".dat");

  /*-- print header --*/

  if (myid == 0) {
    thefile = fopen (filename, "a");
    fprintf(thefile,  "%%\n%% 1.time  2.eta_avg  3.psi_avg  4.N  5.E_kin");
    fprintf(thefile,  "  6.E_pot  7.No   8.A\n%%\n");
    fclose(thefile); 
  }
}


/*-----------------------------------------------------------*/

void info_output(int grid, double t)
{
  FILE          *thefile;

  double         psi, psisq;
  double         No, vortices, Ek, Ep, A;
  double         u, v, sq, max;

  double         sumu     = 0;
  double         sumv     = 0;
  double         mx       = 0;
  double         sumsq    = 0;
  double         sumAk    = 0;
  
  int            i, n, nv, N;

  N = N0*pow(2,grid);
  n = N*N/np;

  /*-- local sums --*/ 

  for (i=0; i<n; i++) {

      u = data[i][0];
      v = data[i][1];
      sq = u*u + v*v;

      sumu    += u;
      sumv    += v;
      sumsq   += sq;
 
      if (sq>mx) mx=sq;


      u = fdata[i][0];
      v = fdata[i][1];
 
      sumAk += u*u + v*v;
      
  }

  mx = sqrt(mx);

  /*-- global sums --*/

  MPI_Reduce(&sumu,    &u,       1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumv,    &v,       1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumsq,   &psisq,   1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumAk,   &A,       1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Allreduce(&mx,   &max,     1,  MPI_DOUBLE, MPI_MAX,    MPI_COMM_WORLD);


  /*-- rescale to represent total or average quantities --*/

  psisq   =  psisq * L*L/(N*N);
  u       =  u / (N*N);
  v       =  v / (N*N);
  No      =  u*u+v*v;

  A       =  A * L*L/(N*N);

  
  /*-- output to file --*/

  Ek = Ekin(grid);
  Ep = 0;

  if (myid == 0) {

    thefile = fopen (filename, "a");
    fprintf(thefile, "%18.12e  %18.12e  %18.12e  %18.12e  ", 
	    t, u, v, psisq);  
    fprintf(thefile, "%18.12e  %18.12e  %18.12e  %18.12e\n", 
	    Ek, Ep, No, A);  
    fclose(thefile); 

  }



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

