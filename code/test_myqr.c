/*
 cc myqr.c  test_myqr.c -o test_myqr.x
*/


#include <stdio.h>


extern void myqr(double *b, double *x);
/* ---------------------------------------------------------- */
int
main(int argc, char **argv)
{
  double f[9];
  double s[5];

  double x0=0.2; 
  double y0=-0.3;
  double c = 1;
  double h=4;

  double dx, dy, r;
  int i, j, k;

  /*-- prepare RHS --*/

  k=0;
  for (j=-1; j<2; j++) {
    for (i=-1; i<2; i++)  {
      dx = i-x0;
      dy = j-y0;
      f[k] = h - 0.5*c*(dx*dx + dy*dy);
      k++;
    }
  }


  /*-- solve linear system --*/

  myqr(f, s);

  printf("\n Solution of linear system and RMS residue:\n\n");
  for (k=0; k<5; k++) printf("%8.4f\n", s[k]);

  /*-- interpretation of results --*/

  printf("\n Collapse parameters:\n\n");

  printf(" Input:  x0 = %f,   y0 = %f,  h = %f,  c = %f\n", x0,y0,h,c);

  c = -2*s[3];
  x0 = s[1]/c;
  y0 = s[2]/c;

  h = s[0] + 0.5*c*(x0*x0 + y0*y0);

  printf(" Output: x0 = %f,   y0 = %f,  h = %f,  c = %f\n\n", x0,y0,h,c);

 
  return 0;

}


/* ----------------------------------------------------- */
