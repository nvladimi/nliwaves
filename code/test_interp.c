#include "header.h"

const  int     N  = 256;
const  double  dx = 0.1;
const  int     gp = 3;


int main(int argc, char **argv)
{
  fftw_complex 	**psi, **psigp, **dd,  **pxxgp, **pyygp, **work;

  ctrl_str  ctrl;
  geom_str  geom;
  phys_str  phys;
  diag_str  diag;
  ic_str    ic;

  collapse_str theC;


  double    L, x, y, h, r;
  int       I,J;

  /* --------------------------------------------- */

  MPI_Init(&argc, &argv);

  /*-- parameters --*/

  L = dx*N;
  x = L/3;
  y = L/6;
  r = 0.1*L;
  h = 1;

  /*-- location of coarse maximum --*/

  theC.i = floor(x/dx + 0.5);
  theC.j = floor(y/dx + 0.5);


  /*-- fill geometry and initial conditions structures --*/

  geom.dx = dx;
  geom.dy = dx;
  geom.Nx = N;
  geom.Ny = N;
  geom.Lx = L;
  geom.Ly = L;
  geom.gp = gp;

  //sprintf(ic.type, "gauss"); 
  strcpy(ic.type, "gauss");

  ic.BumpHeight = h;
  ic.BumpRadius = r;
  ic.BumpX      = x;
  ic.BumpY      = y;
 


  /*-- major array allocation --*/

  arr2D_init(&geom);

  psi     = arr2D_create(0);
  psigp   = arr2D_create(gp);

  dd      = arr2D_create(0);
  pxxgp   = arr2D_create(gp);
  pyygp   = arr2D_create(gp);

  work    = arr2D_create(0);

  /*-- initialization --*/

  fft_init(&geom, work[0]);
  c_interpolate_init(&geom);
  interpol_init(pxxgp, pyygp);

  ic_set(psi, &geom, &ic);



  /*-- prepare second derivatives --*/

  arr2D_copy_gp(psi, psigp);

  arr2D_copy(psi, dd);
  fft_deriv_xx(dd[0]);
  arr2D_copy_gp(dd, pxxgp);

  arr2D_copy(psi, dd);
  fft_deriv_xx(dd[0]);
  arr2D_copy_gp(dd, pyygp);


  /*-- interpolate --*/
 
  c_interpolate(&theC, psigp);

  printf("\n");
  printf("  Input Gaussian:  h=%f, r=%f, x=%f, y=%f\n", h, r, x, y);
  printf("  Second derivative: d2f/dr2 = -2h/rr = %f\n", -2*h/(r*r));

  printf("  Output:\n");
  printf("    x = %f,  y = %f\n", (theC.i+theC.x)*dx, (theC.j+theC.y)*dx);
  printf("    h = %f (least squares),  h = %f (bicubic)\n", theC.h, theC.h0);
  printf("    ddh = %f (least squares),  ddh = %f (bicubic)\n", 
	 theC.ddh, theC.ddh0);
  printf("    ddu = %f,  ddv=%f\n", theC.ddu, theC.ddv);
  printf("\n");
 

  MPI_Finalize();
  exit(0); 
}
