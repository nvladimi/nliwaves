/*
 cc test_dgels.c -o test_dgels.x -llapack -lblas -lg2c
*/


#include <stdio.h>


static long
dgtsv(long N, long NRHS, double *DL, double *D, double *DU, double *B,
      long LDB)
{
  extern void dgtsv_(const long *Np, const long *NRHSp, double *DL,
		     double *D, double *DU, double *B, const long *LDBp,
		     long *INFOp);
  long info;
  dgtsv_(&N, &NRHS, DL, D, DU, B, &LDB, &info);
  return info;
}




static long
dgels(char TRANS, long M, long N, long NRHS, 
      double *A, long LDA, double *B, long LDB, 
      double *WORK, long LWORK)
{
  extern void dgels_(const char *TRANSp, const long *Mp, const long *Np, const long *NRHSp, 
		     double *A, const long *LDAp, double *B, const long *LDBp, 
		     double *WORK, const long *LWORKp, long *INFOp);
  long info;
  dgels_(&TRANS, &M, &N, &NRHS, A, &LDA, B, &LDB, WORK, &LWORK, &info);
  return info;  
}


/* ----------------------------------------------------- */

int
main_simple()
{
  const long n = 3;
  const long m = 4;

  double A[] = {1,0,0,1, 0,2,0,2, 0,0,4,4};
  double f[] = {1,1,1,1};

  double w[50];

  int info, k;

  info = dgels('N', m, n, 1, A, m, f, m, w, 50);
  if (info != 0) fprintf(stderr, "failure with error %d\n", info);

  for (k=0; k<m; ++k) printf("%8.4f\n", f[k]);

  return 0;

}

/* ---------------------------------------------------------- */
int
main(int argc, char **argv)
{
  //a=[1,-1,-1,2;1,0,-1,1;1,1,-1,2;1,-1,0,1;1,0,0,0;1,1,0,1;1,-1,1,2;1,0,1,1;1,1,1,2]

  const long n = 4;
  const long m = 9;

  double A[36];
  double f[9];

  double w[50];

  double x0=0.2; 
  double y0=-0.3;
  double c = 1;
  double h=4;

  double dx, dy;

  int i, j, k, I, J;
  int info;

  k=0;
  for (j=-1; j<2; j++) {
    for (i=-1; i<2; i++)  {

      A[k]   = 1;
      A[k+m] = i;
      A[k+2*m] = j;
      A[k+3*m] = i*i + j*j;

      dx = i-x0;
      dy = j-y0;
      f[k] = h - 0.5*c*(dx*dx + dy*dy);

      printf("%5.1f %5.1f %5.1f %5.1f   %8.4f\n", 
      	     A[k], A[k+n], A[k+2*n], A[k+3*n], f[k] );
      k++;
    }
  }

  info = dgels('N', m, n, 1, A, m, f, m, w, 50);
  if (info != 0) fprintf(stderr, "failure with error %d\n", info);

  // rows 1 to n of B contain the least squares solution vectors; 
  // the residual sum of squares for the solution in each column 
  // is given by the sum of squares of elements N+1 to M in that column; 

  printf("\n Solution:\n");
  for (k=0; k<n; ++k) printf("%8.4f\n", f[k]);

  printf("\n Residual:\n");
  for (k=n; k<m; ++k) printf("%8.4f\n", f[k]);

  return 0;


}



