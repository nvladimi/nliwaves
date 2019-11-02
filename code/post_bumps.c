#include "header.h"

/*---------------------------------------------------------------*/
/*
   The program 

   a) read array of NP, H, and blurred NP
   b) find local maxima from blurred NP
   c) computes H(r) and NP(r) around each maximum
   d) write binary output in the form

      for each maximum, one row contaning
       
       x, y, R_critical,  N(1), N(2), ..., N(nr)  

      or

       x, y, R_critical,  H(1), H(2), ..., H(nr)  

*/

/*----------------------------------------------------------------*/

static int      myid, np;
static int      Nx, Ny;
static float    dx;

static int      nr;    /* max radius in pixels, also GP */

void  radsum(float **A, float **Ar, float Acr, int i0, int j0, int nbump);
void laplace(float **psi, int *In, int *Jn, float *L, int ntot);
void output_txt(float **Nrad, float **Hrad, 
		float *Ln, float *LBn, int ntot, char *filename);

float** alloc_float(int nx, int ny);

/* -------------------------------------------------------------- */

void post_bumps(post_ptr post, geom_ptr geom)

{

  float       **psi,  **ener;  /* field data                      */
  float       **buff;          /* scratch buffer                  */
  float       **Nrad, **Hrad;  /* radial dependence for each bump */
  float        *Ln, *LBn;      /* laplacians at the centers       */
  int          *In, *Jn;       /* index coordiantes of bumps      */

  int           n, nmax, ntot, fn;

  char          filename[80];

  /*-- initialization --*/

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  post_io_init();

  Nx = geom->N0;
  Ny = geom->N0/np;
  dx = geom->Lx/Nx;

  nr   = ceil(post->rbumps / dx); 
  nmax = post->nbumps;

  psi      = alloc_float(Nx+2*nr, Ny+2*nr);
  ener     = alloc_float(Nx+2*nr, Ny+2*nr);
  buff     = alloc_float(Nx+2*nr, Ny+2*nr);
  Nrad     = alloc_float(nr+3, nmax);
  Hrad     = alloc_float(nr+3, nmax);

  In       =  (int *) malloc( nmax * sizeof(int));
  Jn       =  (int *) malloc( nmax * sizeof(int));
  LBn      =  (float *) malloc( nmax * sizeof(float));
  Ln       =  (float *) malloc( nmax * sizeof(float));

  /*-- main loop --*/

  for (fn=post->fn_beg; fn<=post->fn_end; fn+=post->fn_inc) {


    /*-- read blurred psi --*/

    sprintf(filename, "%s.Nblur%04d.%04d", post->runname, 
    	    (int)(post->rblur[0]*1000), fn);

    io_read_float_gp(psi[0],  buff[0], Nx, Ny, nr, filename);


    ntot = find_centers(psi, In, Jn, nmax, post->hbumps);


    /*-- compute lapacian of blurred psi at the centers --*/

    laplace(psi, In, Jn, LBn, ntot);


    /*-- read un-blurred psi and energy --*/

    sprintf(filename, "%s.Nblur%04d.%04d", post->runname, 0, fn);
    io_read_float_gp(psi[0],  buff[0], Nx, Ny, nr, filename);


    sprintf(filename, "%s.Hblur%04d.%04d", post->runname, 0, fn);
    io_read_float_gp(ener[0], buff[0], Nx, Ny, nr, filename);



    /*-- compute lapacian of un-blurred psi at the centers --*/

    laplace(psi, In, Jn, Ln, ntot);


    /*-- compute profiles around each center --*/

    for (n = 0; n < ntot ; n++){

      radsum(psi,   Nrad, 11.68, In[n], Jn[n], n);
      radsum(ener,  Hrad, 0.0,   In[n], Jn[n], n);

    }

    /*-- print profiles to the binary file --*/

    //printf("nr=%d, ntot=%d\n", nr, ntot);

    sprintf(filename, "%s.Nrad.%04d", post->runname, fn);
    io_write_bin(Nrad[0], (nr+3)*ntot*sizeof(float), filename);

    sprintf(filename, "%s.Hrad.%04d", post->runname, fn);
    io_write_bin(Hrad[0], (nr+3)*ntot*sizeof(float), filename);

    /*-- print summary  to the text file --*/

    sprintf(filename, "TXT/%s.rad.%04d.txt", post->runname, fn);
    output_txt(Nrad, Hrad, Ln, LBn, ntot, filename);


  }

}

/*----------------------------------------------------------------*/

int find_centers(float **psi, int *In, int *Jn, int nmax, double hmin){

  int i, j, ntot;
  double hh;

  hh = hmin*hmin;

  ntot = 0;

  for (j=nr; j<Ny+nr; j++) for (i=nr; i<Nx+nr; i++) {

    if ( (psi[j][i] > hh) && 
         (psi[j][i] >= psi[j-1][i]) && (psi[j][i] > psi[j+1][i]) &&
         (psi[j][i] >= psi[j][i-1]) && (psi[j][i] > psi[j][i+1]) ) {

      //printf("n= %2d:  i= %4d,  j= %4d,  h= %6.2f\n", 
      //	     ntot, i-nr, j-nr,  psi[j][i]);
    
      In[ntot] = i-nr;
      Jn[ntot] = j-nr;

      ntot++;

      if (ntot>nmax) return(-1);

    }
  }

  return(ntot);

}
/*----------------------------------------------------------------*/

void laplace(float **psi, int *In, int *Jn, float *L, int ntot){

  int i, j, n;
  double d, overdxdx;

  overdxdx = 1/(dx*dx);

  for (n=0; n<ntot; n++){

    i = In[n]+nr;
    j = Jn[n]+nr;

    d =  sqrt(psi[j-1][i]) + sqrt(psi[j+1][i]) + 
         sqrt(psi[j][i-1]) + sqrt(psi[j][i+1]) - 4*sqrt(psi[j][i]);

    L[n] = d*overdxdx;

  }

}

/*----------------------------------------------------------------*/

void  radsum(float **A, float **Ar, float Acr, int i0, int j0, int nbump){


  float   *Atot;
  int     i,j,k, kk, kkmax;
  float   a, b, rcr, dxdx;

  Ar[nbump][0] = i0*dx;
  Ar[nbump][1] = (myid*Ny+j0)*dx;
 
  Atot = Ar[nbump]+3;

  i0 = i0 + nr;
  j0 = j0 + nr;

  kkmax = nr*nr;
  rcr = -1;

  for (k=0; k<nr; k++) Atot[k] = 0;


  /*-- find amount of stuff in each ring A(r) --*/

  for (j=0; j<Ny+2*nr; j++) for (i=0; i<Nx+2*nr; i++) {

    kk = (i-i0)*(i-i0) + (j-j0)*(j-j0);

    if (kk >= kkmax) continue;
    
    k = floor(sqrt(kk));

    Atot[k] += A[j][i];

  }


  /*-- take into account the size of the cell  --*/

  dxdx = dx*dx;

  for (k=0; k<nr; k++)  Atot[k] =  Atot[k]*dxdx;


  /*-- integrate from the center; find critical radius --*/

  for (k=1; k<nr; k++) {

    Atot[k] += Atot[k-1];

    if ( (rcr < 0) && (Atot[k] - Acr)*(Atot[k-1] - Acr) <= 0) {
      a = Atot[k-1] - Acr;
      b = Atot[k]   - Acr;
      rcr = ( (k-1) + a/(a-b) )*dx;
    }

  }

  if (rcr < 0) {
    if (Atot[nr-1]-Acr > 0) rcr = 0;
    else rcr = nr*dx;
  }

  Ar[nbump][2] = rcr;

}

/*----------------------------------------------------------------*/

void output_txt(float **Nrad, float **Hrad, 
		float *Ln, float *LBn, int ntot, char *filename) {

  FILE *thefile;
  MPI_Status    status;
  int n, nout;

  /*-- get current output number --*/

  if (myid != 0 ) 
    MPI_Recv(&nout, 1, MPI_INTEGER, myid-1, myid, MPI_COMM_WORLD, &status);
  else {
    thefile = fopen(filename, "wt");
    fprintf(thefile, 
	    "#_1.n  2.x  3.y  4.h  5.rN  6.rH  7.delsq  8.delsq_blurred\n\n");
    fclose(thefile);

    nout = 0;

  }

  /*-- output data --*/


  thefile = fopen(filename, "at");

  for (n=0; n<ntot; n++){
    nout++;
    fprintf(thefile, "%4d  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %12.6e  %12.6e\n",
	    nout, Nrad[n][0], Nrad[n][1], 
	    sqrt(Nrad[n][3])/dx, 
	    Nrad[n][2], Hrad[n][2], 
	    Ln[n], LBn[n]);
  }

  fclose(thefile);


  /*-- communicate current output number --*/

  if (myid != np-1)
    MPI_Send(&nout, 1, MPI_INTEGER, myid+1, myid+1, MPI_COMM_WORLD);

}

/*----------------------------------------------------------------*/

float** alloc_float(int nx, int ny)
{
  int            i, j;
  long int       n, ntot;
  float          **A;

  ntot = nx*ny;

  A = 	  (float **)malloc( ny   * sizeof(float *) );
  A[0] =  (float *) malloc( ntot * sizeof(float));

  for (j=1; j<ny; j++)  A[j] = A[0] + j*nx;

  for (n=0; n<ntot; n++) A[0][n]=0;

  return(A);
}

/*----------------------------------------------------------------*/
