#include "header.h"

static  int   myid, np;

static  int   tmrRK;

extern void   dealias(int grid);
extern void   dealias_init(geom_ptr geom, fftw_complex *psi, fftw_complex *psihat);

extern void   rk_one_step(int grid, double dt);          // Uo := Uo + dU
extern void   rk_init(geom_ptr geom, fftw_complex *psi, int np_in);
extern void   rhs_init(geom_ptr geom, phys_ptr phys, int np_in);

/* ---------------------------------------------------------------- */

void
evolve_init(geom_ptr geom, phys_ptr phys, 
	       fftw_complex *psi, fftw_complex *psihat){

  
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

 
  /*--- create suppementary arrays ---*/

  rk_init(geom, psi, np);
  rhs_init(geom, phys, np); 

  dealias_init(geom, psi, psihat);

  
  /*--- timers ---*/

  tmrRK    = timer_set("RK4");
 
}

/* ---------------------------------------------------------------- */

void evolve_one_step(int grid, double dt)   // Uo := Uo + dU
{

  timer_on(tmrRK);

  rk_one_step(grid, dt);
  
  dealias(grid);

  timer_off(tmrRK);

}

/* ---------------------------------------------------------------- */

void evolve_biglin(int grid,  double dt){
}


/* ---------------------------------------------------------------- */

