#include "header.h"

static  int   tmrRK;

extern void   dealias(int grid);
extern void   dealias_init(geom_ptr geom, fftw_complex *psi, fftw_complex *psihat);

extern void   rk_one_step(int grid, double dt);          // Uo := Uo + dU
extern void   rk_init(geom_ptr geom, fftw_complex *psi);

/* ---------------------------------------------------------------- */

void evolve_init(geom_ptr geom, phys_ptr phys, fftw_complex *psi, fftw_complex *psihat)
{

  rk_init(geom, psi);

  rhs_init(geom, phys, psi); 

  dealias_init(geom, psi, psihat);

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

