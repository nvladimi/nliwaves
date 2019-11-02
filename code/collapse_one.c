#include "header.h"

static  int 		 myid, np;
static  int 		 N0, Nx, Ny;
static  double           Lx, dx;

static  char             clpsfile[80];
static  char             msg[120];

static  double         **C;

static  fftw_complex    *Psi;          // 2D array
static  fftw_complex    *Phi;          // 2D array
static  fftw_complex    *psi1D;        // Slices in X/Y direction
static  fftw_complex    *tmp1D;        // Slices in X/Y direction

const static int      records1 =  6;   //  [h, phase, ddh, ddp, halfwidth, phi]
const static int      records2 = 10;   //  [t, h, phase, phi, 2(ddh, ddp, halfwidth)]

static  int              nsize;
static  int              step;
static  int              current_grid;

static  double           nu;

void process_slice(fftw_complex *psi, double *center);


/*----------------------------------------------------------------*/

int collapses_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag, phys_ptr phys){

  int       Ngrids, n;  
  FILE      *thefile;
  char      line[80];

  if (diag->clpsCmax <= 0) return(0);

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  Ngrids = geom -> Ngrids;
  N0     = geom -> N0;
  Lx     = geom -> Lx;

  nu     = phys ->coefNu;

  nsize  = diag -> clpsNmax;

  strcpy (clpsfile, ctrl->runname);
  strcat (clpsfile, ".clps");


  /*-- allocate work arrays of psi and |psi| with ghost points --*/

  Nx  =  N0 * pow(2, Ngrids-1);
  Ny  =  Nx/np;

  psi1D      = (fftw_complex *) malloc(  Nx * sizeof(fftw_complex));
  tmp1D      = (fftw_complex *) malloc(  Nx * sizeof(fftw_complex));

  Phi        = (fftw_complex *) malloc(  Nx*Ny * sizeof(fftw_complex) );

  C =     (double **)malloc( nsize * sizeof(double *) );
  C[0] =  (double *) malloc( nsize*(records2) * sizeof(double));

  for (n=1; n<nsize; n++)  C[n] = C[0] + n * (records2);

  /*-- print header to the output file --*/

  step = 0;

  if (myid == 0) {
    thefile = fopen (clpsfile, "a");
    fprintf(thefile, 
	    "%%_1.t   2.|psi|  3.phase  4.phi_x  5.6.|psi|_rr  7.8.phase_rr  9.10.halfwidth\n");
    fclose(thefile);
  }


  sprintf(msg, "#msg: collapses are initialized\n");
  io_message(msg,0);

  return(1);

}

/*----------------------------------------------------------------*/

void collapses_prepare(fftw_complex *psiData, int grid)
{

  int i;

  Psi = psiData;
  current_grid = grid;

  for (i=0; i<Nx*Ny; i++) {
    Phi[i][0] = Psi[i][0]*Psi[i][0] + Psi[i][1]*Psi[i][1];
    Phi[i][1] = 0;
  }

  fft_davey_stewartson_phi_x(Phi, grid, nu);

  Nx  =  N0 * pow(2, grid);
  Ny  =  Nx/np;

  dx  = Lx/Nx;

}

/*----------------------------------------------------------------*/
/* Scan the whole domain for local maxima;                        */
/* match found maxima against all collapses in database;          */
/* if match not found, start new history.                         */

void collapses_update(double t)
{
  double        centerX[records1];
  double        centerY[records1];

  int             I, J, i, j, i0;

  MPI_Status    status;
  MPI_Request   request;


  I   = Nx/2;
  J   = Nx/2;

  if (step!=0) if (C[step-1][0] == t) return;

  /*-- x-direction --*/

  if (myid == np/2) {
    if (np == 1)  i0 = Nx*Nx/2; 
    else i0 = 0; 
    for (i=0; i<Nx; i++)  {
      psi1D[i][0] = Psi[i0+i][0];
      psi1D[i][1] = Psi[i0+i][1];
    }
    process_slice(psi1D, centerX);
    centerX[5] = Phi[i0+Nx/2][0];
    if (np>1) 
      MPI_Send(centerX, records1, MPI_DOUBLE, 0, 12, MPI_COMM_WORLD);
  }

  if ((myid == 0) && (np>1))  
  MPI_Recv(centerX, records1, MPI_DOUBLE, np/2, 12, MPI_COMM_WORLD, &status);
  

  /*-- y-direction --*/

  memset(tmp1D, 0, 2*Nx*sizeof(double));

  for (j=0; j<Ny; j++)  {
    tmp1D[myid*Ny + j][0] = Psi[Nx*j + I][0];
    tmp1D[myid*Ny + j][1] = Psi[Nx*j + I][1];
  } 

  MPI_Reduce(tmp1D, psi1D, 2*Nx,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (myid == 0)  process_slice(psi1D, centerY);

 
  /*-- store collapse info in database --*/

  C[step][0] = t;
  C[step][1] = centerX[0];  // h 
  C[step][2] = centerX[1];  // phase
  C[step][3] = centerX[5];  // phi
    
  if (myid == 0) {
    for (i=2; i<records1-1; i++) {
      C[step][2*i]   = centerX[i];
      C[step][2*i+1] = centerY[i];
    }
  }
  

  step++;

}

/*-----------------------------------------------------------*/
/* Print to text file the history of the collapse.           */

void collapses_remove(double t)
{
  int       s, r;
  FILE      *thefile;


  /*-- write to file --*/

  if (myid == 0) {
 
    thefile = fopen (clpsfile, "a");

    for (s=0; s<step; s++){
      for (r=0; r<records2; r++) fprintf(thefile, "  %23.16e", C[s][r]);
      fprintf(thefile, "\n"); 
    }

    fclose(thefile); 
  }

  /*-- clear database --*/

  memset(C[0], 0, step*(records2)*sizeof(double) );

  step = 0;

}

/*-----------------------------------------------------------*/

void collapses_addnew(double t) {}

/*-----------------------------------------------------------*/

void process_slice(fftw_complex *psi, double *center){

  double  u, v, h, x, h1, h2;
  double  ux,  vx,  uxx,  vxx; 
  double  ddu, ddv, ddp, ddh;

  int     I, i;

  I = Nx/2;

  /*-- collapse height --*/

  u = psi[I][0];
  v = psi[I][1];
  h = sqrt(u*u + v*v);

  /*-- half width --*/

  i=I;  h1=h; h2=h;

  while (h2 > 0.5*h) {
    i++;
    h1 = h2;   
    h2 = sqrt(psi[i][0] * psi[i][0]  +  psi[i][1] * psi[i][1]); 
  }

  x = (i - (h2-0.5*h)/(h2-h1) - I) * dx;


  /*-- second derivatives at the center --*/

  fft_deriv_xx_1D(psi, current_grid);

  ddu = psi[I][0];
  ddv = psi[I][1];

  ddp  = (ddv*u - ddu*v)/(h*h);
  ddh  = (ddu*u + ddv*v)/h;


  /*-- save data --*/

  center[0] = h;
  center[1] = atan2(u,v);
  center[2] = ddh;
  center[3] = ddp;
  center[4] = x;

}

/*-----------------------------------------------------------*/
