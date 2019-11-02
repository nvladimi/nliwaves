/*----------------------------------------------------------------*/
/* 
   Miscelanious functions, alternative to those in the code, 
   saved just in case.  Some have issues.

   Do not link in.

*/
/*----------------------------------------------------------------*/


#include "header.h"

#define MM    49
#define MMMM  2401 

static  int 		 myid, np;
static  int 		 Nx, Ny, gp, IP;
static  double           dx;

static  char             msg[80];
static  int              tmrInterpol; 

static  double           Qt[MMMM], Ri[16];

static fftw_complex     **px, **py, **pxx, **pyy, **pxy, **dd;

/*----------------------------------------------------------------*/

extern void deriv_cntdiff( collapse_str *theC, fftw_complex **psi, 
			   int I, int J, double x0, double y0);
extern void qr(double *f, double *s);

extern double bicubic(double *A, double x0, double y0);

extern void c_int_bicubic(collapse_str *theC, fftw_complex **psi);
extern int c_int_bicubic_spectral(collapse_str *theC, fftw_complex **psi);

/*----------------------------------------------------------------*/

void c_interpolate_init(geom_ptr geom){

  int       m, c, i, j, tot;  
  FILE      *thefile;
  char      line[80];

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  Nx =  geom->Nx;
  Ny =  geom->Ny;
  dx =  geom->dx;
  gp =  geom->gp;
  IP = 3;

  /*-- allocating arrays with derivatives --*/

  dd    = arr2D_create(0);
  px    = arr2D_create(gp);
  py    = arr2D_create(gp);
  pxx   = arr2D_create(gp);
  pyy   = arr2D_create(gp);
  pxy   = arr2D_create(gp);


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
/* Prepare spectral derivatives of psi, needed for interpolation */

void c_prepare_deriv(fftw_complex **psi){

  arr2D_copy(psi, dd);
  fft_deriv_x(dd[0]);
  arr2D_copy_gp(dd, px);

  arr2D_copy(psi, dd);
  fft_deriv_y(dd[0]);
  arr2D_copy_gp(dd, py);

  arr2D_copy(psi, dd);
  fft_deriv_xx(dd[0]);
  arr2D_copy_gp(dd, pxx);

  arr2D_copy(psi, dd);
  fft_deriv_yy(dd[0]);
  arr2D_copy_gp(dd, pyy);

  arr2D_copy(psi, dd);
  fft_deriv_xy(dd[0]);
  arr2D_copy_gp(dd, pxy);

}

/*----------------------------------------------------------------*/
/* Interpolate properties of the current collapse from field data */

void c_interpolate( collapse_str *theC, fftw_complex **psi)
{
  int adj;

  timer_on(tmrInterpol);

  adj=1;

  while (adj){

    adj=c_int_bicubic_spectral(theC, psi);
  }


  timer_off(tmrInterpol);
}

/*----------------------------------------------------------------*/
/*
   The center of collapse is found as the point of zero gradient,
   computed with Newton's method.  The derivatives at the grid
   points are computed sectrally.  The radial derivative is 
   approximated as
 
      u_rr = u_xx + u_yy - (x u_x + y u_2 ) / r^2. 

   The final quantities are computed on the grid and interpolated 
   to the center using bicubic interpolation. 

   Note that members of the collapse structure, theC has different 
   meaning now:   0 < x,y < 1 (not -0.5 <x,y <0.5)  and index i,j 
   corresponts to the point to the left of the centre (not closest).

*/

int c_int_bicubic_spectral(collapse_str *theC, fftw_complex **psi){

  const  int    m       = 4;
  const  double tol     = 1.e-4;
  const  int    itermax = 50;

  double  Au[16], Av[16], Ah[16], Addu[16], Addv[16], Addp[16], Addh[16];
  double  DX[16], DY[16];
  double  x0, y0, x, y, rr; 

  double  u, v, h;
  double  ux, uy, vx, vy, uxx, uyy, vxx, vyy, uxy, vxy;
  double  fx, fy, fxx, fyy, fxy, D, dX, dY; 
  double  ddu, ddv, ddp, ddh;
  double  err, errmax;

  int     i, j, i0, j0, I, J, iter, adj;

  I  = theC->i;
  J  = theC->j;


  x0 = theC->x;
  y0 = theC->y;
  i0 = theC->i + gp-1;
  j0 = theC->j + gp-1;



  /*-- prepare to iterate x0,y0 --*/

  for (j=0; j<m; j++) for (i=0; i<m; i++) {

    u = psi[j0+j][i0+i][0];
    v = psi[j0+j][i0+i][1];

    ux  = px[j0+j][i0+i][0];
    vx  = px[j0+j][i0+i][1];
    uy  = py[j0+j][i0+i][0];
    vy  = py[j0+j][i0+i][1];
    uxx = pxx[j0+j][i0+i][0];
    vxx = pxx[j0+j][i0+i][1];
    uyy = pyy[j0+j][i0+i][0];
    vyy = pyy[j0+j][i0+i][1];
    uxy = pxy[j0+j][i0+i][0];
    vxy = pxy[j0+j][i0+i][1];
 
    fx  = u*ux + v*vx;
    fy  = u*uy + v*vy;
    fxx = u*uxx + ux*ux + v*vxx + vx*vx;
    fyy = u*uyy + uy*uy + v*vyy + vy*vy;
    fxy = u*uxy + ux*uy + v*vxy + vx*vy;

    D   = fxx*fyy - fxy*fxy;

    DX[j*m+i] = (fx*fyy - fy*fxy)/D;
    DY[j*m+i] = (fy*fxx - fx*fxy)/D;

  }


  /*-- iterate x0,y0 using Newton's method--*/

  printf(" x0=%f, y0=%f   => \n", x0,y0);

  err = 1.0;
  iter = 0;
  adj = 0;
  errmax=0.1;


  while ((err>tol) && (iter<itermax)) {

    dX = bicubic(DX,   x0, y0)/dx;
    dY = bicubic(DY,   x0, y0)/dx;

    err = sqrt(dX*dX + dY*dY);

    if (err>errmax) {
      dX=dX*errmax; 
      dY=dY*errmax;
    }

    x0 = x0 - dX;
    y0 = y0 - dY;

    iter++;

    printf(" x0=%f, y0=%f,  err=%e, iter=%d\n", x0,y0,err,iter);

    if (x0<0)  {theC->x=x0+1; theC->i--; adj++;}
    if (x0>=1) {theC->x=x0-1; theC->i++; adj++;}
    if (y0<0)  {theC->y=y0+1; theC->j--; adj++;}
    if (y0>=1) {theC->y=y0-1; theC->j++; adj++;}

    if (adj) return(1);

  }

  printf("done\n");
 
  theC->x=x0;
  theC->y=y0;


  /*-- collect quantities to interpolate  --*/

  for (j=0; j<m; j++) for (i=0; i<m; i++) {

    x  = (i-1-x0)*dx;
    y  = (j-1-y0)*dx;

    rr = x*x + y*y;

    u = psi[j0+j][i0+i][0];
    v = psi[j0+j][i0+i][1];
    h = sqrt(u*u + v*v);

    ux  = px[j0+j][i0+i][0];
    vx  = px[j0+j][i0+i][1];
    uy  = py[j0+j][i0+i][0];
    vy  = py[j0+j][i0+i][1];
    uxx = pxx[j0+j][i0+i][0];
    vxx = pxx[j0+j][i0+i][1];
    uyy = pyy[j0+j][i0+i][0];
    vyy = pyy[j0+j][i0+i][1];
    uxy = pxy[j0+j][i0+i][0];
    vxy = pxy[j0+j][i0+i][1];
 
 
    ddu  = uxx + uyy - (x*ux + y*uy)/rr;
    ddv  = vxx + vyy - (x*vx + y*vy)/rr;

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

  /*-- interpolate everything --*/

  theC->u   = bicubic(Au,   x0, y0);
  theC->v   = bicubic(Av,   x0, y0);
  theC->h   = bicubic(Ah,   x0, y0);
  theC->ddu = bicubic(Addu, x0, y0);
  theC->ddv = bicubic(Addv, x0, y0);
  theC->ddp = bicubic(Addp, x0, y0);
  theC->ddh = bicubic(Addh, x0, y0);

  return(0);

}

/*----------------------------------------------------------------*/

void deriv_cntdiff( collapse_str *theC, fftw_complex **psi, 
		    int I, int J, double x0, double y0){

  fftw_complex tmp[4], pxx, pyy;

  double x1,y1, a, b, c1, c2, c3, c4;

  double over36  = 1./36;
  double over6hh = 1./(6*dx*dx);

  int     i, j, i1, j1;


  if (x0>=0) {x1=x0; i1=I+gp-1;} else {x1=1+x0; i1=I+gp-2;}
  if (y0>=0) {y1=y0; j1=J+gp-1;} else {y1=1+y0; j1=J+gp-2;}


  /*-- x-derivatives --*/

  a=x1; b=1-a;
  c1=b;  c2=1-3*b;  c3=1-3*a;  c4=a;

  for (j=0; j<4; j++){

    tmp[j][0] =  c1*psi[j1+j][ i1 ][0] + c2*psi[j1+j][i1+1][0] +
                 c3*psi[j1+j][i1+2][0] + c4*psi[j1+j][i1+3][0];
    tmp[j][1] =  c1*psi[j1+j][ i1 ][1] + c2*psi[j1+j][i1+1][1] +
                 c3*psi[j1+j][i1+2][1] + c4*psi[j1+j][i1+3][1];
  }

  a=y1; b=1-a;
  c1=-a*b*(b+1);  c2=3*b*(a+1)*(b+1);  c2=3*a*(a+1)*(b+1);  c4=-a*b*(a+1);
 
  pxx[0] = (c1*tmp[0][0] + c2*tmp[1][0] + 
	     c3*tmp[2][0] + c4*tmp[3][0])*over6hh;
  pxx[1] = (c1*tmp[0][1] + c2*tmp[1][1] + 
	     c3*tmp[2][1] + c4*tmp[3][1])*over6hh;


  /*-- y-derivatives --*/

  a=y1; b=1-a;
  c1=b;  c2=1-3*b;  c3=1-3*a;  c4=a;

  for (i=0; i<4; i++){

    tmp[i][0] =  c1*psi[ j1 ][i1+i][0] + c2*psi[j1+1][i1+i][0] +
                 c3*psi[j1+2][i1+i][0] + c4*psi[j1+3][i1+i][0];
    tmp[i][1] =  c1*psi[ j1 ][i1+i][1] + c2*psi[j1+1][i1+i][1] +
                 c3*psi[j1+2][i1+i][1] + c4*psi[j1+3][i1+i][1];
  }

  a=x1; b=1-a;
  c1=-a*b*(b+1);  c2=3*b*(a+1)*(b+1);  c2=3*a*(a+1)*(b+1);  c4=-a*b*(a+1);
 
  pyy[0] = (c1*tmp[0][0] + c2*tmp[1][0] + 
	     c3*tmp[2][0] + c4*tmp[3][0])*over6hh;
  pyy[1] = (c1*tmp[0][1] + c2*tmp[1][1] + 
	     c3*tmp[2][1] + c4*tmp[3][1])*over6hh;


  /*-- values --*/

  a=x1; b=1-a;
  c1=-a*b*(b+1);  c2=3*b*(a+1)*(b+1);  c2=3*a*(a+1)*(b+1);  c4=-a*b*(a+1);

  for (j=0; j<4; j++){

    tmp[j][0] =  c1*psi[j1+j][ i1 ][0] + c2*psi[j1+j][i1+1][0] +
                 c3*psi[j1+j][i1+2][0] + c4*psi[j1+j][i1+3][0];
    tmp[j][1] =  c1*psi[j1+j][ i1 ][1] + c2*psi[j1+j][i1+1][1] +
                 c3*psi[j1+j][i1+2][1] + c4*psi[j1+j][i1+3][1];
  }

  a=y1; b=1-a;
  c1=-a*b*(b+1);  c2=3*b*(a+1)*(b+1);  c2=3*a*(a+1)*(b+1);  c4=-a*b*(a+1);
 
  theC->u = (c1*tmp[0][0] + c2*tmp[1][0] + 
            c3*tmp[2][0] + c4*tmp[3][0])*over36;
  theC->v = (c1*tmp[0][1] + c2*tmp[1][1] + 
            c3*tmp[2][1] + c4*tmp[3][1])*over36;

  /*-- update --*/

  theC->ddu = 0.5*(pxx[0] + pyy[0]);
  theC->ddv = 0.5*(pxx[1] + pyy[1]);

}

/*----------------------------------------------------------------*/
/*


    
    ux =  (  psi[j0+j][i0+i-2][0] - 8*psi[j0+j][i0+i-1][0] + 
	   8*psi[j0+j][i0+i+1][0] -   psi[j0+j][i0+i+2][0] )*over12dx;

    uy =  (  psi[j0+j-2][i0+i][0] - 8*psi[j0+j-1][i0+i][0] +
	   8*psi[j0+j+1][i0+i][0] -   psi[j0+j+2][i0+i][0] )*over12dx;

    vx =  (  psi[j0+j][i0+i-2][1] - 8*psi[j0+j][i0+i-1][1] + 
	   8*psi[j0+j][i0+i+1][1] -   psi[j0+j][i0+i+2][1] )*over12dx;

    vy =  (  psi[j0+j-2][i0+i][1] - 8*psi[j0+j-1][i0+i][1] +
	   8*psi[j0+j+1][i0+i][1] -   psi[j0+j+2][i0+i][1] )*over12dx;



    uxx = ( -    psi[j0+j][i0+i-2][0]
	    + 16*psi[j0+j][i0+i-1][0]
	    - 30*psi[j0+j][i0+i  ][0]
	    + 16*psi[j0+j][i0+i+1][0]
	    -    psi[j0+j][i0+i+2][0] )*over12dx2;

    uyy = ( -    psi[j0+j-2][i0+i][0]
	    + 16*psi[j0+j-1][i0+i][0]
	    - 30*psi[j0+j  ][i0+i][0]
	    + 16*psi[j0+j+1][i0+i][0]
	    -    psi[j0+j+2][i0+i][0] )*over12dx2;

    vxx = ( -    psi[j0+j][i0+i-2][1]
	    + 16*psi[j0+j][i0+i-1][1]
	    - 30*psi[j0+j][i0+i  ][1]
	    + 16*psi[j0+j][i0+i+1][1]
	    -    psi[j0+j][i0+i+2][1] )*over12dx2;

    vyy = ( -    psi[j0+j-2][i0+i][1]
	    + 16*psi[j0+j-1][i0+i][1]
	    - 30*psi[j0+j  ][i0+i][1]
	    + 16*psi[j0+j+1][i0+i][1]
	    -    psi[j0+j+2][i0+i][1] )*over12dx2;

    ddu  = uxx + uyy - (x*ux + y*uy)/rr;
    ddv  = vxx + vyy - (x*vx + y*vy)/rr;

*/

/*----------------------------------------------------------------*/

