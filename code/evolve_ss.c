#include "header.h"

static  int 		 myid, np;
static  int 		 N0;
static  double           L;
static  double           pi;

static  double           dealiasZ;
static  int              filtersize;
static  double           *filter;
static  fftw_complex    **force;
static  fftw_complex    **force_add;

static  fftw_complex    *data, *fdata;

static  double           coefA, coefB, coefC, coefE, expoS;
static  int              focus;

static  double           f_kmin, f_kmax, f_fmax;
static  char             f_type[80];

static  int              tmrLin, tmrNonlin;
static  char             msg[80];

static void step_lin(int grid, double dt, int save_psihat);
static void step_nonlin(int grid, double dt);
static void init_force_const(int grid, double a);
static void init_force_DF(double a, double b, double c, double f_kdamp);
static void init_force_hiK(double a, double c);
static void init_force_cascade(double beta, double gamma, double f_kdamp, double f_kfric);

static void update_force(int grid, double dt);
static void update_force_loK(int grid);
static void update_force_add(double dt);


/* ---------------------------------------------------------------- */

void
evolve_init(geom_ptr geom, phys_ptr phys, 
	       fftw_complex *psi, fftw_complex *psihat){

  int     i, j, Nx, Ny, filtersize, Ngrids;
  double  x;

  data = psi;
  fdata = psihat;

  //FILE *thefile;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  N0 =  geom->N0;
  L  =  geom->Lx;
  pi =  acos(0)*2;


  coefA = phys->coefA;
  coefB = phys->coefB;
  coefC = phys->coefC;
  coefE = phys->coefE;
  expoS = phys->expoS;
  focus = phys->focus;


  strcpy(f_type,  phys->f_type);
  f_kmin     =  phys->f_kmin;
  f_kmax     =  phys->f_kmax;
  f_fmax     =  0;


  dealiasZ   =  geom->dealiasZ;
  filtersize =  geom->dealiasF;
  Ngrids     =  geom->Ngrids;

  /*--- set up dealising filter ---*/

  filter = (double *)malloc( filtersize*sizeof(double) );

  for (i=0; i<filtersize; i++){
    x = 1 - 1.*i/filtersize;
    filter[i] = x*x*(3-2*x);
  }

  /*--- create arrays of forcing ---*/

  Nx = N0*pow(2, Ngrids-1);
  Ny = Nx/np;

  force    =  (fftw_complex **)malloc( Ny * sizeof(fftw_complex *) );
  force[0] =  (fftw_complex *) malloc( Nx*Ny * sizeof(fftw_complex));

  for (j=1; j<Ny; j++)  force[j] = force[0] + j*Nx;

  force_add    =  (fftw_complex **)malloc( Ny * sizeof(fftw_complex *) );
  force_add[0] =  (fftw_complex *) malloc( Nx*Ny * sizeof(fftw_complex));

  for (j=1; j<Ny; j++)  force_add[j] = force_add[0] + j*Nx;

  memset(force[0],     0, Nx*Ny*sizeof(fftw_complex));
  memset(force_add[0], 0, Nx*Ny*sizeof(fftw_complex));


  /*--- can use force[j][i] only for single grid setups ---*/


  if( strcmp(phys->f_type, "none") == 0 ) {
    init_force_const(Ngrids-1, 0);
  }

  else if( strcmp(phys->f_type, "mult_const") == 0 ){
    init_force_const(Ngrids-1, coefB*coefE);
  }

  else if( strcmp(phys->f_type, "mult_DF") == 0 ) {
    if (Ngrids != 1) {
      sprintf(msg, "# ERROR: Forcing works only on a single grid. Exiting.\n");
      io_message(msg,0);
      exit(-1);
    }
    init_force_DF(phys->f_alpha, phys->f_beta, phys->coefB, phys->f_kdamp);
  }


  else if( strcmp(phys->f_type, "add_DF") == 0 ) {
    if (Ngrids != 1) {
      sprintf(msg, "# ERROR: Forcing works only on a single grid. Exiting.\n");
      io_message(msg,0);
      exit(-1);
    }
    init_force_DF(0, phys->f_beta,  0, phys->f_kdamp);
    f_fmax = phys->f_alpha;
    srand48(myid);
  }


  else if( strcmp(phys->f_type, "add_cascade") == 0 ) {
    if (Ngrids != 1) {
      sprintf(msg, "# ERROR: Forcing works only on a single grid. Exiting.\n");
      io_message(msg,0);
      exit(-1);
    }
    init_force_cascade(phys->f_beta, phys->coefB, phys->f_kdamp, phys->f_kfric);
    f_fmax = phys->f_alpha;
    srand48(myid);
  }



  else if( strcmp(phys->f_type, "mult_hiK") == 0 ) {
   if (Ngrids != 1) {
      sprintf(msg, "# ERROR: Forcing works only on a single grid. Exiting.\n");
      io_message(msg,0);
      exit(-1);
    }
    init_force_hiK(phys->f_alpha, phys->coefB);
  }

  else if( strcmp(phys->f_type, "mult_random") == 0 ){
    f_fmax = phys->f_alpha;
  }

  else {
    sprintf(msg, "# ERROR: unknown forcing %s\n", phys->f_type);
    io_message(msg,0);
    exit(-1);
  }


  //thefile = fopen ("force.tmp", "wb");
  //fwrite (force[0], Nx*Ny*sizeof(double), 1, thefile);
  //fclose(thefile); 


  /*--- timers ---*/

  tmrLin    = timer_set("splitstep:lin");
  tmrNonlin = timer_set("splitstep:non");

}

/* ---------------------------------------------------------------- */

void update_force(int grid, double dt){

  if( strcmp(f_type, "add_DF") == 0 )  update_force_add(dt);
  if( strcmp(f_type, "add_cascade") == 0 )  update_force_add(dt);

}


/* ---------------------------------------------------------------- */

void init_force_const(int grid, double a){

  int i, ntot, N;

  N = N0*pow(2,grid);

  ntot = N*N/np;

  for (i=0; i<ntot; i++) {
    force[0][i][0] = a;
    force[0][i][1] = 0;
  }
}

/* ---------------------------------------------------------------- */

void init_force_DF(double a, double b, double c, double f_kdamp){

  int     i, j, n, kx, ky, N;
  double  pump, damp, q, qq, h, kk;
  double  kkr, kkl;

  N = N0;  /* this forcing does not make sense for adjustable grids */

  n = N/np;

  kkl = f_kmin * f_kmin;
  kkr = f_kmax * f_kmax;


  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    force[j][i][1] = 0;

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    kk = kx*kx + ky*ky;

    if (kk<kkl) {
      force[j][i][0] = c;
      continue;
    }


    if ((kk > kkl) && (kk < kkr))
      pump = sqrt((kk - kkl)*(kkr - kk));
    else
      pump = 0;

    q  = sqrt(kk)/f_kdamp;
    qq = q*q;


    if       (q<1.e-6) h = 0;
    else if  (q<1)     h = exp(5. - 5./qq) /(6*qq*q*qq);
    else               h = 1 - 5/6. * exp(0.5 - 0.5*qq); 

    damp = kk*h;

    force[j][i][0] = a*pump - b*damp;

   // --- HACK!!! Checking out imaginary forcing --- 
   //
   //  force[j][i][0] = - b*damp;
   //  force[j][i][1] =   a*pump;

  }

}

/* ---------------------------------------------------------------- */

static void init_force_cascade(double beta, double gamma, double k_damp, double k_fric) {

  int     i, j, n, kx, ky, N;
  double  kk_damp, kk_fric, kk, k, q, c;


  N = N0;  /* this forcing does not make sense for adjustable grids */
  c = 1;  //c=pow(L/N, 2);  /* empirical coefficient, to be verified */

  n = N/np;

  kk_damp = k_damp * k_damp;
  kk_fric = k_fric * k_fric;

  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    kk = kx*kx + ky*ky;


    if (kk > kk_damp) {

      q  = sqrt(kk)/k_damp;

      q = - beta * c * q*q*q*q* (q-1)*(q-1);

      /* q = - beta * kk *(1 - 5/6. * exp(0.5 - 0.5*q*q));  //  DF forcing */

      force[j][i][0] = q;
      force[j][i][1] = 0;

    }

    if (kk <= kk_fric) {
       force[j][i][0] = -gamma / sqrt(kk);
       force[j][i][1] =  0;
    }

    if (kk == 0) {
       force[j][i][0] = -gamma;
       force[j][i][1] =  0;
    }

  }

}

/* ---------------------------------------------------------------- */

void init_force_hiK(double a, double c){

  int     i, j, n, kx, ky, N;
  double  kk_min, kk, k, q;


  N = N0;  /* this forcing does not make sense for adjustable grids */

  n = N/np;

  kk_min = f_kmin * f_kmin;

  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    force[j][i][1] = 0;

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    kk = kx*kx + ky*ky;

    if (kk <= kk_min) {
      force[j][i][0] = c;
    } else {
      k = sqrt(kk);
      q = (k-f_kmin)/(f_kmax-f_kmin);            /*  0<q<1  */
      force[j][i][0] =  4*a*q*(1 - q);
      /* force[j][i][0] =  4*a*q*q*(1 - q*q); */
      /* alternatively: force[j][i] =  a[1 - (2q-1)^4]; */
    }

  }

}
/* ---------------------------------------------------------------- */

void update_force_loK(int grid){

  int     i, j, n, kx, ky, N;
  double  kk_min, kk_max, kk, k, q, ph;


  N = N0*pow(2, grid);
  n = N/np;

  kk_min = f_kmin * f_kmin;
  kk_max = f_kmax * f_kmax;

  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    kk = kx*kx + ky*ky;

    if ( (kk>kk_min) && (kk<kk_max) ) {
      k = sqrt(kk);
      q = (k-f_kmin)/(f_kmax-f_kmin);            /*  0<q<1  */
      q = f_fmax*q*(1 - q);
      ph = pi*(drand48() - 0.5);
      force[j][i][0] =  q*cos(ph);
      force[j][i][1] =  q*sin(ph);
    } else {
      force[j][i][0] = 0;
      force[j][i][1] = 0;
    }
  }

}

/* ---------------------------------------------------------------- */

void update_force_add(double dt){

  int     i, j, n, kx, ky, N;
  double  kkr, kkl, kk, dkk;
  double  ph, amp, alf, c;

  N = N0;            /* this forcing does not make sense for adjustable grids */
  c = N*N*0.11;     /* semi-empirical cooefficient, might be a wrong guess */
  // c = N*N/pi/pi;


  n = N/np;

  kkl = f_kmin * f_kmin;
  kkr = f_kmax * f_kmax;

  dkk = kkr - kkl;

  alf = c * 6 * f_fmax * dt / ( pi * dkk*dkk*dkk);      

  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    kk = kx*kx + ky*ky;

    if ((kk > kkl) && (kk < kkr)) {
      amp = sqrt( alf * (kk - kkl)*(kkr - kk) );
      ph = 2*pi*drand48();
      force_add[j][i][0] =  amp*cos(ph);
      force_add[j][i][1] =  amp*sin(ph);
    }

  }

}

/* ---------------------------------------------------------------- */

void evolve_one_step_ss2nd(int grid, double dt){

  update_force(grid, dt);

  step_lin(grid, dt/2, 0);
  step_nonlin(grid, dt);
  step_lin(grid, dt/2, 1);

}

/* ---------------------------------------------------------------- */

void evolve_one_step(int grid, double dt){

  const double  d1 = 1 / ( 2 - pow(2, 1/3.) );
  const double  d3 = d1;
  const double  d2 = 1 - (d1+d3);

  const double  c1 = 0.5 * d1;
  const double  c2 = 0.5 - c1;
  const double  c3 = 0.5 - c1;
  const double  c4 = c1;

  update_force(grid, dt);

  step_nonlin(grid, c1*dt);
  step_lin   (grid, d1*dt, 0);
  step_nonlin(grid, c2*dt);
  step_lin   (grid, d2*dt, 0);
  step_nonlin(grid, c3*dt);
  step_lin   (grid, d3*dt, 1);
  step_nonlin(grid, c4*dt);

}

/* ---------------------------------------------------------------- */

void step_lin(int grid,  double dt, int save_psihat){

  double         c, d, DD, coefEA;
  double         reE, imE, reQ, imQ, sqQ;
  double         re, im, rePsi, imPsi, reF, imF;
  int            N, n, i, j, kx, ky, dkx, dky;
  int            kmin, kmax;

  timer_on(tmrLin);

  N = N0 * pow(2,grid);

  if ((dealiasZ <= 0) || (dealiasZ >= 1))  kmin = N/2;
  else                                     kmin = N/2*dealiasZ;
  kmax = kmin + filtersize;
  if (kmax > N/2) kmax = N/2; 


  n = N/np;
  c = 4*pi*pi/(L*L);

  coefEA = coefE*coefA;

  fft_wrap(data, grid, FORWARD);
 
  for (j=0; j<n; j++) for (i=0; i<N; i++) {

    kx = i;
    ky = j + myid*n;

    if (kx > N/2) kx = kx-N;
    if (ky > N/2) ky = ky-N;

    if ((abs(kx) < kmax) && (abs(ky) < kmax)) {

      DD =  -(kx*kx + ky*ky)*c;

      rePsi = data[N*j+i][0];
      imPsi = data[N*j+i][1];


      /*--- Q = ((a*epsilon + i)*DD  + forcing )*dt  ---*/

      reQ = (DD*coefEA + force[j][i][0])*dt;
      imQ = (DD        + force[j][i][1])*dt;

      sqQ = reQ*reQ + imQ*imQ; 

      if (sqQ == 0) {
	rePsi +=  force_add[j][i][0]*dt;
	imPsi +=  force_add[j][i][1]*dt;
        goto dealiasing;
      }


      /*--- E = exp(Q) ---*/
      
      d   =   exp(reQ);
      reE = d*cos(imQ);
      imE = d*sin(imQ);


      /*--- F = force_add / Q ---*/

      re  =  force_add[j][i][0];
      im  =  force_add[j][i][1];

      reF = (re*reQ + im*imQ)/sqQ;
      imF = (im*reQ - re*imQ)/sqQ;


      /*--- psi = psi * E ---*/

      re = rePsi;
      im = imPsi;

      rePsi = re*reE - im*imE;
      imPsi = re*imE + im*reE;


      /*--- psi = psi + F*(E-1) ---*/

      rePsi += reF*(reE-1) - imF*imE;
      imPsi += reF*imE     + imF*(reE-1);


    dealiasing:

      /*--- smoothing dealiasing ---*/

      dkx = abs(kx) - kmin;
      dky = abs(ky) - kmin;

      if (dkx >= 0) {
	rePsi = rePsi * filter[dkx];
	imPsi = imPsi * filter[dkx];
      }

      if (dky >= 0) {
	rePsi = rePsi * filter[dky];
	imPsi = imPsi * filter[dky];
      }

    } else {
      rePsi = 0;
      imPsi = 0;
    }

    data[N*j+i][0] = rePsi;
    data[N*j+i][1] = imPsi;

  }

  /*-- at last substep, save a copy of psihat to work array --*/

  if (save_psihat)  fft_save_fourier(data, fdata, grid);


  fft_wrap(data, grid, BACKWARD);

  timer_off(tmrLin);

}
/* ---------------------------------------------------------------- */

void step_nonlin(int grid, double dt){

  fftw_complex   Q, F;
  double         u, v, pp;
  double         a, c, q, ps;
  double         over2EC, overSEC, halfS, s1, s2;
  long int       i, ntot;
  int            N;

  timer_on(tmrNonlin);

  N = N0 * pow(2,grid);

  ntot = N*N/np;

  a = focus*(2+expoS)*coefE*coefC*dt;

  over2EC = 1/(2*coefE*coefC);
  overSEC = 1/(expoS*coefE*coefC);
  s1      = 1/(expoS+2);
  s2      = expoS/(expoS+2);
  halfS   = expoS/2;

  /*---  loop over data points ---*/

  for (i=0; i<ntot; i++) {
 
    u = data[i][0];
    v = data[i][1]; 
 
    pp = u*u + v*v;

    /*--- compute Q ---*/

    if (coefC == 0) {

      Q[0] = 0;
      Q[1] = focus*pp*dt; 

    } else if (coefC < 0){      /* hack for nonlinear saturation, s=2 */

      Q[0] = 0;
      Q[1] = focus*pp*dt* ( 1 + coefC * pp );
 
    } else if (expoS == 0){

      q    =  log(1 + a*pp);
      Q[0] = -q/2;
      Q[1] =  q*over2EC;


    } else {

      if (pp == 0) {Q[0] = 0; Q[1] = 0;}
      else{
	ps   = pow(pp, halfS);                   /*  ps = |psi|^s  */    
	q    = log(1 + a*pp*ps);
	Q[0] = -q*s1;                            /*  s1 = 1/(s+2)  */
	Q[1] =  (exp(q*s2) - 1)*overSEC/ps;      /*  s2 = s/(s+2)  */
      }

    }

    /*--- F = exp(Q) ---*/

    c    =   exp(Q[0]);
    F[0] = c*cos(Q[1]);
    F[1] = c*sin(Q[1]); 


    /*--- psi = psi*F ---*/

    data[i][0] = u*F[0] - v*F[1];
    data[i][1] = u*F[1] + v*F[0];

  }
  
  timer_off(tmrNonlin);

}

/* ---------------------------------------------------------------- */

void evolve_biglin(int grid,  double dt){

  step_lin(grid, dt, 1);

}
/* ---------------------------------------------------------------- */
