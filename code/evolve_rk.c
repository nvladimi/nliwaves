#include "header.h"

static  int   tmrRK;

/* ---------------------------------------------------------------- */

void evolve_init(geom_ptr geom, phys_ptr phys, fftw_complex *psi, fftw_complex *psihat)
{

  rk_init(geom, psi);

  rhs_init(geom, phys, psi); 

  forcing_init(geom, phys, psi, psihat);

  tmrRK    = timer_set("RK4");
 
}

/* ---------------------------------------------------------------- */

void evolve_one_step(int grid, double dt)   // Uo := Uo + dU
{

  timer_on(tmrRK);

  rk_one_step(grid, dt);
  
  forcing(grid, dt);

  timer_off(tmrRK);

}

/* ---------------------------------------------------------------- */

void evolve_biglin(int grid,  double dt){
}


/* ---------------------------------------------------------------- */

