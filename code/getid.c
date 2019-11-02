
/*---------------------------------------------------------------*/
/*                                                               */
/*                                                               */
/* Program to determine process ID for submitting several serial */
/* jobs with PBS script.  See "pequena_serial.sub" for example.  */
/*                                                               */
/*  Compiling on Pequena:                                        */
/*                                                               */
/*    module load mvapich-intel                                  */
/*    mpicc getid.c -o getid.x                                   */
/*                                                               */
/*                                                               */
/*---------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  int id;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  printf("%d", id);
  MPI_Finalize();
  exit(0);
}
