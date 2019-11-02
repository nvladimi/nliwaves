#include "header.h"

static const   int DEBUG = 0;
static const   int gp    = 4;


static const   int EMPTY     = -1;
static const   int ACTIVE    = -2;
static const   int FINISHED  = -3;
static const   int LOST      = -4;


static  int 		 myid, np;
static  int 		 N0, Nx, Ny;
static  double           Lx, dx;
static  double           psi0;
static  double           rate;

static  int              csize;   // array size for collapses stored
static  int              nsize;   // array size for history
static  int              cmax;    // current maximum number of collapses
static  int              nmax;    // current maximum length of history
static  int              nout;    // current output

static  fftw_complex   **psi;       // psi without ghost points
static  fftw_complex   **psigp;     // psi with ghost points
static  double         **Psigp;     // absolute value of psi
static  fftw_complex    *psigpData; // storage for psipg
static  double          *PsigpData; // storage for Psigp

static  collapse_str   **C;       // array of collapse histories, C[c][n]
static  int             *cn;      // array with numbers of records, cn[c]
static  int             *cstatus; // array collapse status, cstatus[c]
static  double          *chmax;   // collapse max heights, chmax[c]
static  double          *cxmax;   // collapse max x-position, cxmax[c]
static  double          *cymax;   // collapse max y-position, cymax[c]
static  double          *ctmax;   // times of collapses are at max, ctmax[c]


static  int              tmrAddnew, tmrUpdate, tmrRemove, tmrPrepare; // timers

static  char             clpsfile[80];
static  char             debugfile[80];
static  char             msg[120];
static  FILE             *dbg;

/*----------------------------------------------------------------*/

static void cdb_record(collapse_str theC, int c, double t);
static int  cdb_match(double x, double y);
static int  cdb_adjust(collapse_str *theC);
static void cdb_output(int c);

/*----------------------------------------------------------------*/

int collapses_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag, phys_ptr phys){

  int       c, j, tot, Ngrids;  
  FILE      *thefile;
  char      line[80];

  if (diag->clpsCmax <= 0) return(0);

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);


  Ngrids = geom -> Ngrids;
  N0     = geom -> N0;
  Lx     = geom -> Lx;

  csize  = diag -> clpsCmax;
  nsize  = diag -> clpsNmax;
  psi0   = diag -> clpsPsi0;
  rate   = diag -> clpsRate;

  strcpy (clpsfile, ctrl->runname);
  strcat (clpsfile, ".clps");

  if (DEBUG) sprintf(debugfile,"debug.%04d", myid);
  

  /*-- allocate array of collapse histories, C[c][n]  --*/

  C = 	  (collapse_str **)malloc( csize * sizeof(collapse_str *) );
  C[0] =  (collapse_str *) malloc( csize*nsize * sizeof(collapse_str));

  for (c=1; c<csize; c++)  C[c] = C[0] + c*nsize;


  /*-- allocate collapse database --*/

  cstatus = (int *)malloc( csize * sizeof(int) );
  cn      = (int *)malloc( csize * sizeof(int) );
  chmax   = (double *)malloc( csize * sizeof(double) );
  cxmax   = (double *)malloc( csize * sizeof(double) );
  cymax   = (double *)malloc( csize * sizeof(double) );
  ctmax   = (double *)malloc( csize * sizeof(double) );

  for (c=0; c<csize; c++)  {
    cstatus[c] = EMPTY;
    cn[c]      = -1;
    chmax[c]   =  0;
    ctmax[c]   =  0;
  }


  /*-- allocate work arrays of psi and |psi| with ghost points --*/

  Nx  =  N0 * pow(2, Ngrids-1);
  Ny  =  Nx/np;
  tot = (Nx+2*gp)*(Ny+2*gp);

  psigpData = (fftw_complex *) malloc( tot * sizeof(fftw_complex));
  PsigpData = (double *) malloc( tot * sizeof(double));

  psi       = (fftw_complex **) malloc(  Ny * sizeof(fftw_complex *));
  psigp     = (fftw_complex **) malloc( (Ny+2*gp) * sizeof(fftw_complex *));
  Psigp     = (double **) malloc( (Ny+2*gp) * sizeof(double *));

  /*-- set the collapse counter to zero --*/

  cmax  = 1;
  nmax  = 0;
  nout  = 0;

  /*-- print header to the output file --*/

  if (myid == 0) {
    thefile = fopen (clpsfile, "a");
    fprintf(thefile, "%%_1.t   2.|psi|  3.|psi|_rr  4.phase_rr  5.phase");
    //fprintf(thefile, "  5.I  6.J  7.x  8.y  9.h0  10.ddh0  11.err");
    //fprintf(thefile, "  12.u  13.v  14.ddu  15.ddv\n");
    fclose(thefile);
  }

  /*-- init array handling (ghostpoint copy) --*/

  arr2D_init(geom, gp);

  /*-- init interpolation part --*/

  cInt_init();


  /*-- timers --*/

  tmrAddnew  = timer_set("clps:addnew");
  tmrUpdate  = timer_set("clps:update");
  tmrRemove  = timer_set("clps:remove");
  tmrPrepare = timer_set("clps:prepare");

  sprintf(msg, "#msg: collapses are initialized\n");
  io_message(msg,0);

  return(1);

}

/*----------------------------------------------------------------*/
/* Both "addnew" and "update" need psi with ghost points.         */
/* If we want interpolate spectral derivatives, compute them here.*/
/* Currently we compute derivatives with finite differences.      */

void collapses_prepare(fftw_complex *psiData, int grid)
{

  int           i,j;
  int           n, buffsize;
  MPI_Status    status;
  MPI_Request   request;

  timer_on(tmrPrepare);

  Nx  =  N0 * pow(2, grid);
  Ny  =  Nx/np;

  dx  = Lx/Nx;

  for (j=0; j<Ny;      j++)  psi[j]   = psiData   + j* Nx;
  for (j=0; j<Ny+2*gp; j++)  psigp[j] = psigpData + j*(Nx+2*gp);
  for (j=0; j<Ny+2*gp; j++)  Psigp[j] = PsigpData + j*(Nx+2*gp);


  arr2D_copy_gp( psi, psigp, grid);
 //misc_prepare_deriv(psi);


  cInt_prepare(psigp, Nx, Ny, gp, dx);

  timer_off(tmrPrepare);

}

/*----------------------------------------------------------------*/
/* Scan the whole domain for local maxima;                        */
/* match found maxima against all collapses in database;          */
/* if match not found, start new history.                         */

void collapses_addnew(double t)
{
  collapse_str     theC;
  int              i, j, tot, c;
  double           x, y, u, v, h;

  timer_on(tmrAddnew);

  /*-- compute Psigp = |psigp| --*/

  tot = (Nx+2*gp)*(Ny+2*gp);

  for (i=0; i<tot; i++) {
    u = psigp[0][i][0]; 
    v = psigp[0][i][1];
    Psigp[0][i] = sqrt(u*u + v*v);
  }

  /*-- find local maxima above threshold and process them --*/

  for (j=gp; j<=Ny+gp; j++)  for (i=gp; i<=Nx+gp; i++) {

    h = Psigp[j][i];

    if ( ( h > psi0) && 
	 ( h >= Psigp[j-1][i]) && ( h > Psigp[j+1][i]) &&
	 ( h >= Psigp[j][i-1]) && ( h > Psigp[j][i+1]) ) {

      x = (i-gp)*dx;
      y = (j-gp + myid*Ny)*dx;

      c = cdb_match(x, y);

 
      /*-- if collapse has no history, add it --*/

      if (cstatus[c] == EMPTY) {

        theC.x = x;
        theC.y = y;
        theC.h = h;

	if (DEBUG) {
	  dbg=fopen(debugfile, "a");
	  fprintf(dbg, "\nADDING     c=%d       (%8.4f, %8.4f)  at t=%f, h=%f\n", 
		  c, theC.x, theC.y, t, theC.h); 
	  fclose(dbg);
	}

	if (cInt_update(&theC) > 0) cdb_record(theC, c, t);

      }

    }
  }

  timer_off(tmrAddnew);

}


/*----------------------------------------------------------------*/
/* Interpolate all active collapses.  Do not attempt to find      */
/* new ones, unless somebody shifted out of bounds.               */

void collapses_update(double t){

  collapse_str     theC;
  int              c, n, flag, err, errall;
  double           dt, h;

  timer_on(tmrUpdate);

  /*-- adjust collapse grid locations --*/

  err = 0;

  for (c=0; c<cmax; c++){

    if (cstatus[c] != ACTIVE) continue;

    n = cn[c];

    if (C[c][n].t == t) continue;
 
    h  = C[c][n].h;
    dt = t -  C[c][n].t;

    if (dt < 1./(h*h*rate))   continue;

    theC.x = C[c][n].x;
    theC.y = C[c][n].y;


    if (DEBUG) {
      dbg=fopen(debugfile, "a");
      fprintf(dbg, "updating   c=%d, n=%d (%8.4f, %8.4f)  at t=%f, h=%f\n", 
	      c, n, theC.x, theC.y, t, h); 
      fclose(dbg);
    }

    flag = cInt_update(&theC);

    if (flag >  0) cdb_record(theC, c, t);  // updated successfully
    if (flag == 0) cstatus[c] = LOST;       // lost 
    if (flag <  0) {                        // shifted out of bounds   
      cstatus[c] = LOST;
      err++;
    }

  }


  /*-- if collapse has shifted out of bounds rescan the whole domain  -- */
  
  MPI_Allreduce(&err, &errall, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if (errall != 0) {
    sprintf(msg, "#cDB: unscheduled rescan of collapses at t = %18.12e\n", t);
    collapses_addnew(t);
    io_message(msg,0);
  }
  

  timer_off(tmrUpdate);

}


/*----------------------------------------------------------------*/
/* Make non-active, save to file, and remove from database        */
/* collapses with height below thereshold "phi0".                 */

void collapses_remove(double t)
{
  int           n, c;
  MPI_Status    mpi_status;
  double        h, h0, hmax;

  timer_on(tmrRemove);

  /*-- output non-active collapses and clean database --*/

  if (myid != 0 ) 
    MPI_Recv(&nout, 1, MPI_INTEGER, myid-1, myid, MPI_COMM_WORLD, &mpi_status);

  for (c=0; c<cmax; c++)  {

    if (cstatus[c] == EMPTY)  continue;
    if (cstatus[c] == ACTIVE) continue;

    /*-- finished or lost --*/

    n    = cn[c];
    h    = C[c][n].h;
    h0   = C[c][0].h;
    hmax = chmax[c];

    /* output only long lived collapses with maximum */
    if ((n >= 8) && (h<hmax) && (h0<hmax)) {
      cdb_output(c); 
      nout++;
    }

    if (DEBUG) {
      dbg=fopen(debugfile, "a");
      fprintf(dbg,"removing   c=%d, n=%d (%8.4f, %8.4f)  at t=%f, h=%f\n", 
	      c, n, C[c][n].x, C[c][n].y, t, C[c][n].h);
      fclose(dbg);
    }

    cstatus[c] = EMPTY; 
    ctmax[c]   = 0; 
    chmax[c]   = 0; 
    cn[c]      = -1;
  }

  /*-- communicate current output number --*/

  if (myid != np-1)
    MPI_Send(&nout, 1, MPI_INTEGER, myid+1, myid+1, MPI_COMM_WORLD);
  else
    MPI_Send(&nout, 1, MPI_INTEGER, 0,      0,      MPI_COMM_WORLD);

  if (myid == 0)
    MPI_Recv(&nout, 1, MPI_INTEGER, np-1,   0, MPI_COMM_WORLD, &mpi_status);

  timer_off(tmrRemove);
}

/*----------------------------------------------------------------*/
/* Record properties of the current collapse in the database.     */

void cdb_record(collapse_str theC, int c, double t)
{
  int n;

  n = cn[c];
  if (n < nsize-1){
    n++;
    cn[c] = n;
  }  
  else {
    sprintf(msg, "#err(%d): collapse history is too short\n",myid);
    io_message(msg,1);
  }

  theC.t = t;

  if (theC.h > chmax[c]){
    chmax[c] = theC.h;
    cxmax[c] = theC.x;
    cymax[c] = theC.y;
    ctmax[c] = t;
  }

  if (( theC.h < 0.5*chmax[c]) || (theC.h < psi0) ) 
    cstatus[c] = FINISHED; 
  else
    cstatus[c] = ACTIVE;


  C[c][n] = theC;
  if (n>nmax) nmax=n;

  if (DEBUG) {
    dbg=fopen(debugfile, "a");
    fprintf(dbg,"recorded   c=%d, n=%d (%8.4f, %8.4f)  at t=%f\n", 
  	   c, cn[c], C[c][n].x, C[c][n].y, t); 
    fclose(dbg);
  }

}

/*----------------------------------------------------------------*/
/* Find collapse in the database which matches given location.    */

int cdb_match(double x, double y)
{
  int    c, n, s;
  double d;

  /*-- find collapse that matched current --*/

  for (c=0; c<cmax; c++)  {

    if (cstatus[c] == ACTIVE){

      n = cn[c];

      d = abs(C[c][n].x - x) + abs(C[c][n].y - y);

      if (d <= 2*dx) return(c);

    }
  }


  /*-- if no match found, use empty slot --*/

  for (c=0; c<csize; c++)  {

    if (cstatus[c] == EMPTY) {

      if (c==cmax) cmax=c+1;    
      return(c);
    
    }
  }

  /*-- if we reached this point, no slots available --*/ 

  sprintf(msg, "#err(%d): number of active collapses exceeded\n",myid);
  io_message(msg,1);

  return(csize-1);

}


/*-----------------------------------------------------------*/
/* Print to text file the history of the collapse.           */

void cdb_output(int c)
{
  int       n, ntot;
  FILE      *thefile;

  ntot = cn[c];
 
  sprintf(msg,
	  "#cDB(%4d):  collapse %4d   cmax = %6d   nmax = %6d   t = %18.12e\n",
	  myid, nout, cmax, nmax, C[c][ntot-1].t);
  io_message(msg,1);

  /*-- write to file --*/

  thefile = fopen (clpsfile, "a");

  fprintf(thefile, "\n\n %%collapse  %4d:", nout);
  fprintf(thefile, "   h_max( %18.12e )= %18.12e", ctmax[c], chmax[c]);
  fprintf(thefile, "   x,y= %10.4e %10.4e", cxmax[c], cymax[c]);
  fprintf(thefile, "   h_end( %10.4e )= %10.4e\n", 
	  C[c][ntot-1].t, C[c][ntot-1].h);



  for (n=0; n<ntot; n++){

    fprintf(thefile, " %19.12e %19.12e %19.12e %19.12e %19.12e\n",
	    C[c][n].t,  C[c][n].h, C[c][n].ddh, C[c][n].ddp,
	    atan2(C[c][n].v, C[c][n].u) );
 

    /*   
    fprintf(thefile, "  %4d  %4d  %6.3f  %6.3f",
	    C[c][n].i,  C[c][n].j + myid*Ny, C[c][n].x,  C[c][n].y);
    fprintf(thefile, "  %12.6e  %12.6e  %12.6e",
	    C[c][n].h0,  C[c][n].ddh0,  C[c][n].err);
    fprintf(thefile, "  %12.6e  %12.6e  %12.6e  %12.6e\n",
	    C[c][n].u,  C[c][n].v,  C[c][n].ddu,  C[c][n].ddv);
    */

  }

  fclose(thefile); 

}

/*-----------------------------------------------------------*/
