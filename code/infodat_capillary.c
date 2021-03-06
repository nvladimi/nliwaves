#include "header.h"

/* Note: sumA_sq  = (1/2) int( eta^2 + psi^2 ) dxdy                    */
/* Note: sumAk_sq = sum(|Ak|^2) = (2/L^2) sumA_sq                      */
/* Note: Spectrum computes average of |Ak|^2 in rings                  */
/* Note: H = (1/2) int( |\nabla eta|^2 + (1+eta) |nabla psi|^2 ) dxdy  */ 


static  int 		 myid, np;
static  int 		 N0;
static  double           L;
static  double           pi;

static  fftw_complex    *data;
static  fftw_complex    *fdata;

static  char             filename[80];
static  char             msg[80];

static  void momentum(int grid, double *Nk, double *Px, double *Py);

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
    fprintf(thefile,  "%%\n%% 1.time  2.E_pot  3.E_kin  4.E_nl  ");
    fprintf(thefile,  "5.eta_rms  6.eta_min  7.eta_max  8.sumA_sq  ");
    fprintf(thefile,  "9.sumAk_sq  10.Px  11.Py \n%%\n");
    fclose(thefile); 
  }
}


/*-----------------------------------------------------------*/

void info_output(int grid, double t)
{
  FILE          *thefile;

  double         psi, psisq, urms;
  double         Ek, Ep, Enl, A, Ak, px, py;
  double         u, v, max, min;

  double         sumu     = 0;
  double         sumv     = 0;
  double         mx       = 0;
  double         mn       = 0;
  double         sumsq    = 0;
  double         sumusq   = 0;
  double         sumAk    = 0;
  
  int            i, n, nv, N;

  N = N0*pow(2,grid);
  n = N*N/np;

  /*-- local sums --*/ 

  for (i=0; i<n; i++) {

      u =  data[i][0];
      v =  data[i][1];

      sumu    += u;
      sumv    += v;
      sumusq  += u*u;
      sumsq   += u*u + v*v;

      if (u > mx) mx = u;
      if (u < mn) mn = u;
      
  }


  /*-- global sums --*/

  MPI_Reduce(&sumu,    &u,       1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumv,    &v,       1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumsq,   &psisq,   1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumusq,  &urms,    1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Allreduce(&mx,   &max,     1,  MPI_DOUBLE, MPI_MAX,    MPI_COMM_WORLD);
  MPI_Allreduce(&mn,   &min,     1,  MPI_DOUBLE, MPI_MIN,    MPI_COMM_WORLD);


  /*-- rescale to represent total or average quantities --*/

  A    =  psisq / (N*N) * L*L* 0.5;  // =  Ak * L*L/2;
  u    =  u / (N*N);
  v    =  v / (N*N);
  urms =  sqrt(urms / (N*N));

  
  rhs_hamiltonian(grid, &Ep, &Ek, &Enl);

  Ek  = Ek*L*L;
  Ep  = Ep*L*L;
  Enl = Enl*L*L;

  momentum(grid, &Ak, &px, &py);

  px = px * (2*pi/L);
  py = py * (2*pi/L);
  
  if (myid == 0) {

    thefile = fopen (filename, "a");
    fprintf(thefile, "%19.12e %19.12e %19.12e %19.12e ",  t, Ep, Ek, Enl);  
    fprintf(thefile, "%19.12e %19.12e %19.12e %19.12e ", urms, min, max, A);  
    fprintf(thefile, "%19.12e %19.12e %19.12e\n", Ak, px, py);  
    fclose(thefile); 

  }

}

/*-----------------------------------------------------------*/

void momentum(int grid, double *Nk, double *Px, double *Py){

  double        a, b, c, sumnk, sumpx, sumpy, nk, px, py;
  int           N, n, i, j, kx, ky;
 
  N = N0*pow(2,grid);

  n = N/np;

  sumpx = 0;
  sumpy = 0;
  sumnk = 0;
  
  /*-- sum over domain --*/
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    a = fdata[N*j+i][0];
    b = fdata[N*j+i][1];

    nk = a*a + b*b;

    sumpx += nk*kx;
    sumpy += nk*ky;

    sumnk += nk;
    
  }


  /*-- global sums --*/

  MPI_Reduce(&sumnk,   &nk,      1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumpx,   &px,      1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumpy,   &py,      1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  c =  1.0 /(N*N);
  
  *Nk   =  nk * c; 
  *Px   =  px * c; 
  *Py   =  py * c; 

}

/*-----------------------------------------------------------*/
