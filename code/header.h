#include <stdio.h>
#include <stdlib.h>
#include <fftw3-mpi.h>
#include <mpi.h>
#include <math.h>
#include <string.h>

#define FORWARD    1
#define BACKWARD  -1

typedef struct ic_param
{
  char   type[32];
  double BumpHeight;
  double BumpRadius;
  double BumpX;
  double BumpY;
  double BumpA;
  double BumpB;
  int    NumberOfBumps;
  int    Seed;
  double Chirp;
  double kmin, kmax;
  double numParticles;
  double nCondensate;
} ic_str, *ic_ptr;

typedef struct geom_param   // space-related parameters
{
  int    Ngrids;            // max number of grids
  int    grid;              // grid number at init/restart time 
  int    N0;                // number of points on coarsest grid, N0xN0  
  double Lx;                // size of domain in x-direction
  double Ly;                // size of domain in y-direction, Ly=Ly
  double dealiasZ;          // relative K, where dealiasing starts
  int    dealiasF;          // size of dealiasing filter, in points
  double regridZup;         // relative K, where regrid-up cond. is checked 
  double regridZdn;         // relative K, where regrid-down cond. is checked 
  double regridTh;          // regridding (up or down) threshold 
  int    sieve;             // size of coarse mesh saved, sieve x sieve points 
  int    klow;              // size of low k part of spectrum saved
} geom_str, *geom_ptr;


typedef struct ctrl_param   // time-related parameters
{
  char   runname[80];       // base for input/output file names
  int    restart;           // file number to restart from
  int    n0;                // coarsest timestep number at restart 
  double CFL;               // dt = CFL * dx^2
  double dt0;               // override dt, on coarsest grid
  double tmax;              // time at the end of simulation
  int    dnInfo;            // frequency of basic output
  int    dnDiag;            // frequency of computing diagnostics 
  int    dnData;            // frequency of diagnostics and binaty output 
  int    njig;              // frequency of turning off nonlinearity
  int    noff;              // number of steps nonlinearity is off               
} ctrl_str, *ctrl_ptr;


typedef struct diag_param
{
  double pdfMax;
  int    pdfNbins;
  int    clpsNmax;
  int    clpsCmax;
  double clpsPsi0;
  double clpsRate;
  double vortexTh;
} diag_str, *diag_ptr;


typedef struct phys_param
{
  double coefA;
  double coefB;
  double coefC;
  double coefE;
  double coefNu;
  double coefRho;
  double expoS;
  int    focus;
  char   f_type[32];
  double f_kmin;
  double f_kmax;
  double f_kdamp;
  double f_kfric;
  double f_alpha;   // fiber: strength of potential 
  double f_beta;    // fiber: steepness of potential 
  double f_radius;  // fiber: radius
} phys_str, *phys_ptr;


typedef struct
{
  double t;
  double x, y; 
  double h, ddh, ddp;
  double u, ddu;
  double v, ddv;
  double h0, ddh0, err;
} collapse_str;


typedef struct
{
  char   type[80];
  char   runname[80];
  int    fn_beg;
  int    fn_end;
  int    fn_inc;
  int    bin_out;
  int    time_slices;
  double dtime;
  int    nblur;
  double rblur[16];
  int    nbumps;
  double rbumps;
  double hbumps;
} post_str, *post_ptr;





/* arrays.c */

extern void     arr2D_init(geom_ptr geom, int gp_in);

extern fftw_complex** arr2D_create(int grid, int gp);

extern void    arr2D_copy(fftw_complex **A, fftw_complex **B, int grid);
extern void    arr2D_copy_gp(fftw_complex **A, fftw_complex **Agp, int grid);
extern void    arr2D_copy_nogp(fftw_complex **Agp, fftw_complex **A, int grid);
extern void    arr2D_free(fftw_complex **A);
extern void    arr2D_deriv_x(fftw_complex **Agp, fftw_complex **B, int grid);
extern void    arr2D_deriv_y(fftw_complex **Agp, fftw_complex **B, int grid);
extern void    arr2D_deriv_xx(fftw_complex **Agp, fftw_complex **B, int grid);
extern void    arr2D_deriv_yy(fftw_complex **Agp, fftw_complex **B, int grid);


/* parameters.c */

extern void    parameters(int *argc, char ***argv, 
			  geom_ptr geom, ctrl_ptr ctrl, 
			  phys_ptr phys, diag_ptr diag);

/* io.c */

extern void    io_init(ctrl_ptr ctrl, geom_ptr geom);
extern void    io_message(char *msg, int mode);
extern void    io_save_tag(int nio, int n0, int grid, double t);
extern int     io_read_data(fftw_complex *psi, int grid);
extern int     io_save_data(fftw_complex *psi, int grid);
extern void    io_save_sieve(fftw_complex *psi, int grid);
extern void    io_save_klow(fftw_complex *fA, int grid);
extern void    io_save_qflux(double *fA, int nio);

/* ic.c */

extern void    ic_set(fftw_complex *psi, geom_ptr geom,  ic_ptr ic);


/* grids.c */

extern void    grids_init(geom_ptr geom, fftw_complex *work);
extern void    fft_regrid_up (fftw_complex *A, int grid);
extern void    fft_regrid_dn (fftw_complex *A, int grid);


/* fftw2_wrap.c OR fftw3_wrap.c  */

extern void    fft_init      (geom_ptr geom, fftw_complex *work);
extern void    fft_wrap      (fftw_complex *A, int grid, int dir);
extern void    fft_wrap_1D   (fftw_complex *A, int grid, int dir);
extern void    fft_time_init (post_ptr post);
extern void    fft_time_wrap (fftw_complex *A, int dir);


/* fft_extra.c */

extern void    fft_save_fourier (fftw_complex *A, fftw_complex *fA, int grid);
extern void    fft_make_fourier (fftw_complex *A, fftw_complex *fA, int grid);

extern void    fft_extra_init  (geom_ptr geom);
extern void    fft_test        (fftw_complex *A, int grid);
extern void    fft_time_test   (fftw_complex *A);
extern void    fft_deriv_x_1D  (fftw_complex *A, int grid);
extern void    fft_deriv_xx_1D (fftw_complex *A, int grid);
extern void    fft_deriv_x     (fftw_complex *A, int grid);
extern void    fft_deriv_y     (fftw_complex *A, int grid);
extern void    fft_deriv_xx    (fftw_complex *A, int grid);
extern void    fft_deriv_yy    (fftw_complex *A, int grid);
extern void    fft_deriv_xy    (fftw_complex *A, int grid);
extern void    fft_laplacian   (fftw_complex *A, int grid);
extern void    fft_blur        (fftw_complex *A, fftw_complex *B,
			      int grid, double rb);

extern void    fft_davey_stewartson_phi_x(fftw_complex *data, 
					  int grid, double nu);



/* evolve_ss.c  OR  evolve_rk.c */

extern void    evolve_init(geom_ptr geom, phys_ptr phys, 
			      fftw_complex *psi, fftw_complex *psihat);

extern void    evolve_one_step(int grid, double dt);
extern void    evolve_biglin(int grid, double dt);


/* infodat.c */

extern void    info_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag,
			 fftw_complex *psi, fftw_complex *psihat);
extern void    info_output(int grid,  double t);
extern double  abspsi_max(int grid);
extern double  Ekin(int grid);


/* diag_nls.c  OR  diag_dse.c */

extern  void   diag_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag,  phys_ptr phys,
			fftw_complex *psi, fftw_complex *psihat);
extern  void   diag_ic(fftw_complex *psi, fftw_complex *psihat,
		       int grid, int nio, double t);
extern  void   diag_compute(fftw_complex *psi, fftw_complex *psihat, int grid);
extern  void   diag_output(int nio);



/* pdf.c */

extern void    pdf_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag,
			fftw_complex *psi);
extern void    pdf_compute(int grid);
extern void    pdf_output(int filecount);
extern void    pdf_output_clean(int filecount);


/* spectrum.c */

extern void    spectrum_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag,
			     fftw_complex *psi, fftw_complex *psihat);
extern void    spectrum_compute(int grid);
extern void    spectrum_output(int filecount);
extern void    spectrum_output_clean(int filecount);

extern int     spectrum_regrid(int grid, int substep, double t);


/* qflux.c */

extern void    qflux_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag, phys_ptr phys,
			     fftw_complex *psi, fftw_complex *psihat);
extern void    qflux_compute(int grid);
extern void    qflux_output(int filecount);
extern void    qflux_output_clean(int filecount);


/* corrfun.c */

extern void    corrfun_init(geom_ptr geom, char *runname, fftw_complex *acf);
extern void    corrfun_compute(int grid);
extern void    corrfun_output(int filecount);

/* qflux.c */

//extern void qflux_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag);
//extern void qflux_compute(fftw_complex *psi, fftw_complex *psihat, int grid);
//extern void qflux_save();


/* collapsesDB.c */

extern int  collapses_init(ctrl_ptr ctrl, geom_ptr geom, diag_ptr diag, phys_ptr phys);
extern void collapses_prepare(fftw_complex *psi, int grid);
extern void collapses_addnew(double t);
extern void collapses_update(double t);
extern void collapses_remove(double t);


/* collapsesInt.c */

extern void cInt_init();
extern void cInt_prepare(fftw_complex **psigp, 
			 int Nx_in, int Ny_in, int gp_in, double dx);

extern int  cInt_update(collapse_str *theC);


/* timers.c */

extern void timers_init(char *filebase, double simtime);
extern void timers_out(double simtime);
extern void timers_reset(double simtime);

extern int  timer_set(char *name);
extern void timer_reset(int id);
extern void timer_on(int id);
extern void timer_off(int id);


/* post_param.c,  post_io.c */

extern void post_param(int *argc, char ***argv, post_ptr post, geom_ptr geom);
extern void post_io_init();
extern int  io_get_grid(char  *filename, int N0);
extern void io_read_bin(void *buff, int buffsize, char *filename);
extern void io_write_bin(void *buff, int buffsize, char *filename);
extern void post_blur(post_ptr post, geom_ptr geom);
extern void post_bumps(post_ptr post, geom_ptr geom);
extern void post_corrfun(post_ptr post, geom_ptr geom);
extern void post_corrfunx(post_ptr post, geom_ptr geom);
extern void post_condensate(post_ptr post, geom_ptr geom);

extern void io_read_float_gp(float *A, float *buff, 
				int Nx, int Ny, int gp, char *filename);

extern int find_centers(float **psi, int *In, int *Jn, int nmax, double hmin);


/* mask.c */

extern void mask_init(geom_ptr geom, phys_ptr phys, fftw_complex *psi);
extern void mask_apply(int grid, double dt);
