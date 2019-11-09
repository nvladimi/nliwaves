#include "header.h"

static fftw_complex   *SendUp, *SendDown, *RecvUp, *RecvDown;
static int      np, myid, idup, iddown;
static int      gp;

static double   Lx;
static int      N0;

/*----------------------------------------------------------------*/

void arr2D_init(geom_ptr geom, int gp_in)
{
  int Nmax, Ngrids, buffsize;

  N0 = geom->N0;
  Lx = geom->Lx;
  gp = gp_in;

  Nmax = N0*pow(2, geom->Ngrids - 1);

  /* processor ID's and border flags */
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (myid != np-1) { idup   = myid+1; }
  else              { idup   = 0;      }
  if (myid != 0 )   { iddown = myid-1; }
  else              { iddown = np-1;   }

  /* create  buffers for local boundaries exchange */
  buffsize = (Nmax+2*gp)*gp;
  SendUp   = (fftw_complex*)malloc( buffsize * sizeof(fftw_complex));
  SendDown = (fftw_complex*)malloc( buffsize * sizeof(fftw_complex));
  RecvUp   = (fftw_complex*)malloc( buffsize * sizeof(fftw_complex));
  RecvDown = (fftw_complex*)malloc( buffsize * sizeof(fftw_complex));

}


/*----------------------------------------------------------------*/

fftw_complex** arr2D_create(int grid, int gpc)
{
  int            i, j;
  long int       Nx, Ny, n, ntot;
  fftw_complex **A, a;


  Nx = N0*pow(2, grid);
  Ny = Nx/np;

  ntot = (Nx+2*gpc)*(Ny+2*gpc);

  A = 	  (fftw_complex **)malloc( (Ny+2*gpc) * sizeof(fftw_complex *) );
  A[0] =  (fftw_complex *) malloc( ntot * sizeof(fftw_complex));

  for (j=1; j<Ny+2*gpc; j++)  A[j] = A[0] + j*(Nx+2*gpc);

  for (n=0; n<ntot; n++) {
    A[0][n][0]=0;
    A[0][n][1]=0;
  }

  return(A);
}


/*----------------------------------------------------------------*/

void arr2D_free(fftw_complex **A)
{
  free(A[0]);
  free(A);
}

/*----------------------------------------------------------------*/

void arr2D_copy(fftw_complex **A, fftw_complex **B, int grid)
{
  int      Nx, Ny, i,j;

  Nx = N0*pow(2, grid);
  Ny = Nx/np;

  for (j=0; j<Ny; j++)  for (i=0; i<Nx; i++)  {
      B[j][i][0] = A[j][i][0];
      B[j][i][1] = A[j][i][1];
    }
}

/*----------------------------------------------------------------*/

void arr2D_copy_gp(fftw_complex **A, fftw_complex **Agp, int grid)
{
  int           Nx, Ny, i,j;
  int           n, buffsize;
  MPI_Status    status;
  MPI_Request   request;

  Nx = N0*pow(2, grid);
  Ny = Nx/np;

  buffsize = (Nx+2*gp)*gp;


 /* copy middle part */
  for (j=0; j<Ny; j++)  for (i=0; i<Nx; i++)  {
      Agp[j+gp][i+gp][0] = A[j][i][0];
      Agp[j+gp][i+gp][1] = A[j][i][1];
  }

  /* fill buffers */
  n = 0;
  for (j=0; j<gp; j++)
    for (i=0; i<Nx+2*gp; i++) {
      SendDown[n][0] = Agp[gp+j][i][0];
      SendDown[n][1] = Agp[gp+j][i][1];
      SendUp[n][0]   = Agp[Ny+j][i][0];
      SendUp[n][1]   = Agp[Ny+j][i][1];
      n++;
    }

  /* exchange between processors */
  MPI_Irecv(RecvDown, 2*buffsize, MPI_DOUBLE, iddown, 3, MPI_COMM_WORLD, &request);
  MPI_Send(SendUp,    2*buffsize, MPI_DOUBLE, idup,   3, MPI_COMM_WORLD);
  MPI_Wait(&request, &status);
  MPI_Irecv(RecvUp,   2*buffsize, MPI_DOUBLE, idup,   4, MPI_COMM_WORLD, &request);
  MPI_Send(SendDown,  2*buffsize, MPI_DOUBLE, iddown, 4, MPI_COMM_WORLD);
  MPI_Wait(&request, &status);

  /* fill ghost points */
  n = 0;
  for (j=0; j<gp; j++)
    for (i=0; i<Nx+2*gp; i++) {
      Agp[j][i][0]    = RecvDown[n][0];
      Agp[j][i][1]    = RecvDown[n][1];
      Agp[Ny+gp+j][i][0] = RecvUp[n][0];
      Agp[Ny+gp+j][i][1] = RecvUp[n][1];
      n++;
    }

  /* periodic in x direction */

  for (i=0; i<gp; i++) {
    for (j=0; j<Ny+2*gp; j++) {
      Agp[j][i][0] = Agp[j][Nx+i][0];
      Agp[j][i][1] = Agp[j][Nx+i][1];
    }
    for (j=0; j<Ny+2*gp; j++) {
      Agp[j][Nx+gp+i][0] = Agp[j][gp+i][0];
      Agp[j][Nx+gp+i][1] = Agp[j][gp+i][1];
    }

    /*
    for (j=0; j<Ny+2*gp; j++) Agp[j][1] = Agp[j][Nx+1];
    for (j=0; j<Ny+2*gp; j++) Agp[j][0] = Agp[j][Nx];

    for (j=0; j<Ny+2*gp; j++) Agp[j][Nx+2] = Agp[j][2];
    for (j=0; j<Ny+2*gp; j++) Agp[j][Nx+3] = Agp[j][3];
    */
  }

}


/*----------------------------------------------------------------*/

void arr2D_copy_nogp(fftw_complex **Agp, fftw_complex **A, int grid)
{
  int      Nx, Ny, i,j;

  Nx = N0*pow(2, grid);
  Ny = Nx/np;

  for (j=0; j<Ny; j++)  for (i=0; i<Nx; i++)  {
      A[j][i][0] = Agp[j+gp][i+gp][0];
      A[j][i][1] = Agp[j+gp][i+gp][1];
    }

}

/*----------------------------------------------------------------*/

void arr2D_deriv_x(fftw_complex **A, fftw_complex **B, int grid)
{
  int      Nx, Ny, i,j;
  double   c, dx;

  Nx = N0*pow(2, grid);
  Ny = Nx/np;

  dx = Lx/Nx;
  c  = 1/(12*dx);

  for (j=gp; j<Ny+gp; j++)  for (i=gp; i<Nx+gp; i++) {

    B[j-gp][i-gp][0] = 
      (A[j][i-2][0] - 8*A[j][i-1][0] + 8*A[j][i+1][0] - A[j][i+2][0])*c;

    B[j-gp][i-gp][1] = 
      (A[j][i-2][1] - 8*A[j][i-1][1] + 8*A[j][i+1][1] - A[j][i+2][1])*c;

  }

}

/*----------------------------------------------------------------*/

void arr2D_deriv_y(fftw_complex **A, fftw_complex **B, int grid)
{
  int      Nx, Ny, i,j;
  double   c, dx;

  Nx = N0*pow(2, grid);
  Ny = Nx/np;

  dx = Lx/Nx;
  c  = 1/(12*dx);

  for (j=gp; j<Ny+gp; j++)  for (i=gp; i<Nx+gp; i++) {

    B[j-gp][i-gp][0] = 
      (A[j-2][i][0] - 8*A[j-1][i][0] + 8*A[j+1][i][0] - A[j+2][i][0])*c;

    B[j-gp][i-gp][1] = 
      (A[j-2][i][1] - 8*A[j-1][i][1] + 8*A[j+1][i][1] - A[j+2][i][1])*c;

  }

}
/*----------------------------------------------------------------*/

void arr2D_deriv_xx(fftw_complex **A, fftw_complex **B, int grid)
{
  int      Nx, Ny, i,j;
  double   c, dx;

  Nx = N0*pow(2, grid);
  Ny = Nx/np;

  dx = Lx/Nx;
  c  = 1/(12*dx*dx);

  for (j=gp; j<Ny+gp; j++)  for (i=gp; i<Nx+gp; i++) {

    B[j-gp][i-gp][0] = (- A[j][i-2][0] + 16*A[j][i-1][0] - 30*A[j][i][0] 
			+ 16*A[j][i+1][0] - A[j][i+2][0])*c;

    B[j-gp][i-gp][1] = (- A[j][i-2][1] + 16*A[j][i-1][1] - 30*A[j][i][1] 
			+ 16*A[j][i+1][1] - A[j][i+2][1])*c;

  }

}

/*----------------------------------------------------------------*/

void arr2D_deriv_yy(fftw_complex **A, fftw_complex **B, int grid)
{
  int      Nx, Ny, i,j;
  double   c, dx;

  Nx = N0*pow(2, grid);
  Ny = Nx/np;

  dx = Lx/Nx;
  c  = 1/(12*dx*dx);

  for (j=gp; j<Ny+gp; j++)  for (i=gp; i<Nx+gp; i++) {

    B[j-gp][i-gp][0] = (- A[j-2][i][0] + 16*A[j-1][i][0] - 30*A[j][i][0]
			+ 16*A[j+1][i][0] - A[j+2][i][0])*c;

    B[j-gp][i-gp][1] = (- A[j-2][i][1] + 16*A[j-1][i][1] - 30*A[j][i][1]
			+ 16*A[j+1][i][1] - A[j+2][i][1])*c;

  }

}

/*----------------------------------------------------------------*/

void print_ar(double **A, int nx, int ny, int gp)
{
  int     i,j, myid, np;
  int     flag = 0;
  MPI_Status    status;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if (myid==0) flag = 1;
  else  MPI_Recv(&flag, 1, MPI_INTEGER, myid-1, myid, MPI_COMM_WORLD, &status);

  printf("Proc %d:\n", myid);
  for (j=gp; j<gp+ny; j++) {
    for (i=gp; i<gp+nx; i++) printf("%E  ", A[j][i]);
    printf("\n");
  }

  if (myid != np-1)
    MPI_Send(&flag, 1, MPI_INTEGER, myid+1, myid+1, MPI_COMM_WORLD);
}

/*----------------------------------------------------------------*/

