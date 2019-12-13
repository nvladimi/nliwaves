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
    fprintf(thefile,  "%%\n%% 1.time  2.eta_avg  3.psi_avg  4.sumA_sq ");
    fprintf(thefile,  " 5.E_kin  6.E_pot  7.sumAk_sq  8.eta_max \n%%\n");
    fclose(thefile); 
  }
}


/*-----------------------------------------------------------*/

void info_output(int grid, double t)
{
  FILE          *thefile;

  double         psi, psisq;
  double         Ek, Ep, A, Ak;
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

      u =  data[i][0];
      v =  data[i][1];
      sq = u*u + v*v;

      sumu    += u;
      sumv    += v;
      sumsq   += sq;
 
      if (fabs(u)>mx) mx=fabs(u);


      u = fdata[i][0];
      v = fdata[i][1];
 
      sumAk += u*u + v*v;
      
  }


  /*-- global sums --*/

  MPI_Reduce(&sumu,    &u,       1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumv,    &v,       1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumsq,   &psisq,   1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumAk,   &Ak,      1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Allreduce(&mx,   &max,     1,  MPI_DOUBLE, MPI_MAX,    MPI_COMM_WORLD);


  /*-- rescale to represent total or average quantities --*/

  A    =  psisq * L*L/(N*N) * 0.5;  // =  Ak * L*L * 0.5;
  u    =  u / (N*N);
  v    =  v / (N*N);

  Ak   =  Ak/(N*N); 
  
  rhs_hamiltonian(grid, &Ek, &Ep);

  Ek = Ek*L*L;
  Ep = Ep*L*L;
  
  if (myid == 0) {

    thefile = fopen (filename, "a");
    fprintf(thefile, "%19.12e %19.12e %19.12e %19.12e ", 
	    t, u, v, A);  
    fprintf(thefile, "%19.12e %19.12e %19.12e %19.12e\n", 
	    Ek, Ep, Ak, max);  
    fclose(thefile); 

  }

}

/*-----------------------------------------------------------*/
