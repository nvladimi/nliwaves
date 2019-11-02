#include "header.h"

#define MM    49
#define MMMM  2401 

static const   int DEBUG = 0;
static const   int IP    = 3;            // interpolation points

static  fftw_complex   **psi;
static  int 		 myid, np;
static  int              Nx, Ny, gp;
static  double           dx;

static  int              tmrInterpol; 

static  double           Qt[MMMM], Ri[16];

static  char             debugfile[80];
static  char             msg[80];
static  FILE             *dbg;



/*----------------------------------------------------------------*/

static void   cInt_bicubic(collapse_str *theC);

static double bicubic(double *A, double x0, double y0);
static void   qr(double *f, double *s);

/*----------------------------------------------------------------*/

void cInt_init(){

  int       m, c, i, j, tot;  
  FILE      *thefile;
  char      line[80];

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);


  if (DEBUG) sprintf(debugfile,"debug.%04d", myid);

  /*-- read Qt and Ri matrices --*/

  m = (IP*2+1)*(IP*2+1);

  if (IP == 1)      thefile = fopen("qr3x3.DAT", "r");
  else if (IP == 2) thefile = fopen("qr5x5.DAT", "r");
  else if (IP == 3) thefile = fopen("qr7x7.DAT", "r"); 

  if (thefile == NULL) {
    sprintf(msg, "#err:  cannot open QR interpolation file!\n");
    io_message(msg,0);
    exit(-1);
  }

  for (j=0; j<m*m; j++) {
    fgets(line, 80, thefile);
    Qt[j] = atof(line);
  }

  for (j=0; j<16;  j++) {
    fgets(line, 80, thefile);
    Ri[j] = atof(line);
  }

  fclose(thefile);

  /*-- timers --*/

  tmrInterpol = timer_set("clps:interpol");

}

/*----------------------------------------------------------------*/

void cInt_prepare(fftw_complex **psigp, int Nx_in, int Ny_in, 
		  int gp_in, double dx_in){

  psi = psigp;

  Nx = Nx_in;
  Ny = Ny_in;
  gp = gp_in;
  dx = dx_in;

}

/*----------------------------------------------------------------*/
/* Update properties of the current collapse from field data:     */
/* (1) update the location of the center by fitting |psi| with    */
/*     paraboloid,                                                */
/* (2) call cInt_bicubic to interpolate psi and derivatives.      */
/*     returns:  1 is case of success,                            */
/*               0 if collapse is lost, and                       */
/*              -1 if collapse has shifted out of process bounds  */

int cInt_update(collapse_str *theC)
{

  double  f[MM];
  double  s[5];
  double  x, y, u, v;
  double  y_lo, y_hi, L;
  int     i, j, I, J, k, attempt;

  const int max_attempts = 3;

  timer_on(tmrInterpol);

  L    = Nx*dx;
  y_lo = Ny*myid*dx;
  y_hi = Ny*(myid+1)*dx;



  /*-- three attempts to get location of the collapse --*/

  attempt = 0;
  while ( attempt < max_attempts ) {

    x = theC->x;
    y = theC->y;

    I = round(x/dx);
    J = round(y/dx) - myid*Ny;

    if (DEBUG) {
      dbg=fopen(debugfile, "a");
      fprintf(dbg, "   interpolating at (%d, %d) \n", I, J); 
      fclose(dbg);
    }


    /*-- prepare RHS for least square matrix --*/

    k=0;
    for (j=-IP; j<=IP; j++) {
      for (i=-IP; i<=IP; i++)  {

	u = psi[J+gp+j][I+gp+i][0];
	v = psi[J+gp+j][I+gp+i][1];
     
	f[k]  = sqrt(u*u + v*v);
	k++;

      }
    }


    /*-- fit |psi| with with paraboloid to find the center --*/

    if (DEBUG) {
      dbg=fopen(debugfile, "a");
      fprintf(dbg, "   calling QR \n"); 
      fclose(dbg);
    }

    qr(f, s);

    theC->h    =  s[0];
    theC->x    = (s[1] + I)*dx;
    theC->y    = (s[2] + J+Ny*myid)*dx;
    theC->ddh  =  s[3]/(dx*dx);
    theC->err  =  s[4];


    /*-- check if the center shifted out of process bounds */

    if ( theC->y < y_lo || theC->y >= y_hi ) {
      if (DEBUG) {
	dbg=fopen(debugfile, "a");
        fprintf(dbg, "      DB units: %8.4f, %8.4f, h=%f - out of bounds\n", 
	      theC->x, theC->y, theC->h); 

	fclose(dbg);
      }
      timer_off(tmrInterpol);
      return(-1);                            // out of bounds
    }

    if ( theC->x <  0 ) theC->x += L;
    if ( theC->x >= L ) theC->x -= L;

    if (DEBUG) {
      dbg=fopen(debugfile, "a");
      fprintf(dbg, "      DB units: %8.4f, %8.4f, h=%f, %f, %f (%d)\n", 
	      theC->x, theC->y, theC->h, s[1], s[2], attempt); 
      fclose(dbg);
    }


    /*-- check if the center shifted closer to another grid point */

    i = round(s[1]);
    j = round(s[2]);

    if ( (abs(i) > 4) || (abs(j) > 4) ) {
      timer_off(tmrInterpol);
      return(0);                             // collapse is lost
    }

    if ((i==0) && (j==0)) break;             // center found
    else attempt++;

  }

  /*-- all attempts to locate the center are done --*/

  if (attempt == max_attempts) {
    timer_off(tmrInterpol);
    return(0);                                 // collapse is lost
  }

  /*-- if we found the center, interpolate all quantities --*/

  if (DEBUG) {
    dbg=fopen(debugfile, "a");
    fprintf(dbg, "   calling bicubic \n"); 
    fclose(dbg);
  }

  cInt_bicubic(theC);

  if (DEBUG) {
    dbg=fopen(debugfile, "a");
    fprintf(dbg, "   done interpolating\n\n"); 
    fclose(dbg);
  }

  timer_off(tmrInterpol);

  return(1);

}

/*----------------------------------------------------------------*/
/* Interpolate collapse-specific quantities at given location.   */

void cInt_bicubic(collapse_str *theC){

  const  int    m = 4;

  double  Au[16], Av[16], Ah[16];
  double  Addu[16], Addv[16], Addp[16], Addh[16];

  double  x0, y0; // x, y, rr; 
  double  over12dx;
  double  over12dx2;

  double  u, v, h;
  double  ux, uy, vx, vy, uxx, uyy, vxx, vyy; 
  double  ddu, ddv, ddp, ddh;

  int     i0, j0;   /*  first point of interpolation square, without GP  */
  int     i,  j;    /*  running index i,j = 0,1,2,3  */
  int     I,  J;    /*  index in psi array, including  GP  */ 

  over12dx  = 1/(12*dx);
  over12dx2 = 1/(12*dx*dx);


  /*-- convert location to bicubic convention:   0 <= x0,y0 < 1 --*/

  x0 = theC->x;
  y0 = theC->y;

  i0 = floor(x0/dx);
  x0 = x0/dx - i0;
  i0 = i0 - 1;

  j0 = floor(y0/dx);
  y0 = y0/dx - j0;
  j0 = j0 - Ny*myid - 1;

  if (DEBUG) {
    dbg=fopen(debugfile, "a");
    fprintf(dbg, "      BC units: %d+%7.4f, %d+%7.4f\n", i0,x0,j0,y0); 
    fclose(dbg);
  }


  /*-- collect quantities to interpolate  --*/

  for (j=0; j<m; j++) for (i=0; i<m; i++) {

    //x  = (i-1-x0)*dx;
    //y  = (j-1-y0)*dx;

    //rr = x*x + y*y;

    I = i0+i+gp;
    J = j0+j+gp;

    u = psi[J][I][0];
    v = psi[J][I][1];
    h = sqrt(u*u + v*v);


    uxx = ( -    psi[J][I-2][0]
	    + 16*psi[J][I-1][0]
	    - 30*psi[J][I  ][0]
	    + 16*psi[J][I+1][0]
	    -    psi[J][I+2][0] )*over12dx2;

    uyy = ( -    psi[J-2][I][0]
	    + 16*psi[J-1][I][0]
	    - 30*psi[J  ][I][0]
	    + 16*psi[J+1][I][0]
	    -    psi[J+2][I][0] )*over12dx2;

    vxx = ( -    psi[J][I-2][1]
	    + 16*psi[J][I-1][1]
	    - 30*psi[J][I  ][1]
	    + 16*psi[J][I+1][1]
	    -    psi[J][I+2][1] )*over12dx2;

    vyy = ( -    psi[J-2][I][1]
	    + 16*psi[J-1][I][1]
	    - 30*psi[J  ][I][1]
	    + 16*psi[J+1][I][1]
	    -    psi[J+2][I][1] )*over12dx2;


    ddu  = 0.5*(uxx + uyy);
    ddv  = 0.5*(vxx + vyy);

    ddp  = (ddv*u - ddu*v)/(h*h);
    ddh  = (ddu*u + ddv*v)/h;


    Au[j*m+i] = u;
    Av[j*m+i] = v;
    Ah[j*m+i] = h;
    Addu[j*m+i] = ddu;
    Addv[j*m+i] = ddv;
    Addp[j*m+i] = ddp;
    Addh[j*m+i] = ddh;

  }

  /*-- save least square result --*/

  theC->h0   = theC->h;
  theC->ddh0 = theC->ddh;

  /*-- interpolate everything --*/

  theC->u   = bicubic(Au,   x0, y0);
  theC->v   = bicubic(Av,   x0, y0);
  theC->h   = bicubic(Ah,   x0, y0);
  theC->ddu = bicubic(Addu, x0, y0);
  theC->ddv = bicubic(Addv, x0, y0);
  theC->ddp = bicubic(Addp, x0, y0);
  theC->ddh = bicubic(Addh, x0, y0);

}


/*----------------------------------------------------------------*/
/* Generic bicubic interpolation.                                 */

double bicubic(double *A, double x0, double y0){

  const  int  m=4;
  const  double W[] = {  0,  2,  0,  0,
                        -1,  0,  1,  0,
                         2, -5,  4, -1,
                        -1,  3, -3,  1};

  double  B[4], X[4], Y[4], XW[4], YW[4];
  double  p;
  int     i, j;

  /*-- create polynomials t^n --*/

  X[0] = 0.5;
  Y[0] = 0.5;  
  for (j=1; j<m; j++){
    X[j]=X[j-1]*x0;
    Y[j]=Y[j-1]*y0;
  }

  /*-- create vectors XW and YW --*/

  for (j=0; j<m; j++) {
      XW[j] = 0;
      YW[j] = 0;
      for (i=0; i<m; i++) {
	XW[j] += X[i]*W[j+i*m];
	YW[j] += Y[i]*W[j+i*m];
      }	
    }

  /*-- interpolation in x-direction --*/

  for (j=0; j<m; j++) {
    B[j] = 0;
    for (i=0; i<m; i++) B[j] += XW[i]*A[j*m+i];
  }

  /*-- interpolation in y-direction --*/

  p=0;
  for (j=0; j<m; j++) p += YW[j]*B[j];

  return(p);

}

/*-----------------------------------------------------------*/

void
qr(double *b, double *x)
{
  int i,j,m,n; 
  double d[MM];
  double rr, dd, err, c;

  m = (2*IP+1)*(2*IP+1);
  n = 4;

  /* d = Qt*b */
  
  for (j=0; j<m; j++) {
    d[j] = 0;
    for (i=0; i<m; i++) d[j] += Qt[j*m + i]*b[i];
  }


  /* x = Rt*d1 */
  
  for (j=0; j<n; j++) {
    x[j] = 0;
    for (i=0; i<n; i++) x[j] += Ri[j*n + i]*d[i];
  }


  /* residue and relative error */
  rr=0; dd=0;
  for (j=0; j<n; j++) dd +=  d[j]*d[j];
  for (j=n; j<m; j++) rr +=  d[j]*d[j];
  err = sqrt(rr/(rr+dd));


  /* reinterpreting the results */

  c  = 2*x[3];

  x[1] = -x[1]/c;
  x[2] = -x[2]/c;
  x[0] = x[0] - 0.5*c*(x[1]*x[1] + x[2]*x[2]);
  x[3] = c;
  x[4] = err;

}
/*-----------------------------------------------------------*/
