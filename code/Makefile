# Gibbs:     module load openmpi/1.6.5/intel/12.1.5 fftw/3.3.3/openmpi/1.6.5/intel/12.1.5 
# Pequena:   module load openmpi/1.6.3/intel/11.1 fftw/3.3.3/openmpi/1.6.3/intel/11.1
# Encanto:   login to service2; module load mvapich_intel cc/10.1.015

# Name of the executable
NLS     = nls.x
IC      = nls_ic.x
DSE     = nls_ds.x
NLSTEST = nlstest.x
NLSPOST = nlspost.x
CAPILLARY = capillary.x

# List of the sources
# For tests use "test_timers.c", "test_ffttime.c", or "test_deriv.c" as main

C_IC      = main_ic.c grids.c fftw_wrap.c fft_extra.c io.c ic.c  timers.c 

C_NLS     = main.c grids.c fftw_wrap.c fft_extra.c io.c parameters.c evolve_ss.c \
            infodat.c diag_nls.c \
            diag_spectrum.c diag_pdf.c diag_qflux_avg.c\
            timers.c\
            arrays.c  collapseDB.c  collapseInt.c  mask.c

C_DSE     = main.c grids.c fftw_wrap.c  fft_extra.c io.c parameters.c evolve_rk.c \
            rk.c rhs_dse.c dealias_rk.c \
            infodat.c  diag_dse.c diag_spectrum.c \
            timers.c \
            arrays.c  collapse_one.c mask.c

C_NLSPOST = post.c grids.c fftw_wrap.c fft_extra.c corrfun.c \
            post_io.c post_param.c \
            post_blur.c post_bumps.c post_corrfun.c post_corrfunx.c \
            post_condensate.c \
            timers.c

C_NLSTEST = test_deriv.c grids.c fftw_wrap.c fft_extra.c arrays.c io.c ic.c timers.c

C_CAPILLARY = capillary.c grids.c fftw_wrap.c fft_extra.c arrays.c io.c ic.c timers.c


# Compilers
CC   = mpicc
LINK = mpicc
OPT  =  -g3  #-O2


#-------------------------------------------------

#generic

DEF_NLS     = 
#DEF_NLS    = -DDEBUG
#DEF_NLS    = -DQ_FLUX_FILTER
#DEF_NLS    = -DEXIT_AT_PSIMAX



LIB_MPI      = -lmpi
LIB_FFT      = -lfftw3_mpi -lfftw3
INC_MPI      = 
INC_FFT      = 
LIB_ADD      = -lm -lrt



ifeq ($(HOSTNAME), tangle)
        LIB_MPI      = -lmpi
        INC_MPI      =
        INC_FFT      = -I/usr/local/include
        LIB_ADD      = -lm
	DEF_NLS      =
endif



ifeq ($(HOSTNAME), hooloovoo)
        LIB_MPI      = -lmpi
        LIB_FFT      = -L/usr/local/fftw/fftw-3.3.3/lib -lfftw3_mpi -lfftw3_omp -lfftw3 -fopenmp
        INC_MPI      =
        INC_FFT      = -I/usr/local/fftw/fftw-3.3.3/include
        LIB_ADD      = -lm
endif


ifeq ($(HOSTNAME), odin)
        OPT =  -flto -Ofast
endif



ifeq ($(HOSTNAME), mallorn)
        LIB_MPI      = -lmpi
        LIB_FFT      = -L/usr/local/fftw/fftw-2.1.5/lib -lfftw_mpi -lfftw
        INC_MPI      = 
        INC_FFT      = -I/usr/local/fftw/fftw-2.1.5/include
        LIB_ADD      = -lm
endif


ifeq ($(HOSTNAME), athelas)
	LIB_MPI      = -lmpi
	LIB_FFT      = -lfftw_mpi
	INC_MPI      = -I/usr/lib/mpich/include
	INC_FFT      =
endif

ifeq ($(HOSTNAME), cy-c3.lanl.gov)
	LIB_MPI      = -L$(MPIHOME)/lib64 -lmpi 
	LIB_FFT      = -L$(FFTW_HOME)/lib -lfftw_mpi -lfftw
	INC_MPI      = -I$(MPIHOME)/include
	INC_FFT      = -I$(FFTW_INCLUDE)  
endif

ifeq ($(HOSTNAME), pequena)
	LIB_MPI      = 
	LIB_FFT      = -lfftw3_mpi -lfftw3
	INC_MPI      =
	INC_FFT      =
endif


ifeq ($(HOSTNAME), gibbs.alliance.unm.edu)
        LIB_MPI      =
        LIB_FFT      = -lfftw3_mpi -lfftw3
        INC_MPI      =
        INC_FFT      =
endif



ifeq ($(HOSTNAME), stampede)
        LIB_MPI = 
	INC_FFT = -I$(TACC_FFTW3_INC)
	LIB_FFT = -Wl,-rpath,$(TACC_FFTW3_LIB) -L$(TACC_FFTW3_LIB) -lfftw3_mpi -lfftw3
        #LIB_FFT = -L$(TACC_FFTW3_LIB) -lfftw3_mpi -lfftw3
endif


ifeq ($(HOSTNAME), comet-ln2.sdsc.edu)
         #LIB_MPI =
         INC_FFT = -I/opt/fftw/3.3.4/intel/mvapich2_ib/include
         LIB_FFT = -L/opt/fftw/3.3.4/intel/mvapich2_ib/lib  -lfftw3_mpi -lfftw3
endif





#-------------------------------------------------

OBJ_NLS     = $(C_NLS:.c=.o ) $(F_NLS:.f=.o )
OBJ_DSE     = $(C_DSE:.c=.o )
OBJ_IC      = $(C_IC:.c=.o )
OBJ_NLSTEST = $(C_NLSTEST:.c=.o )
OBJ_NLSPOST = $(C_NLSPOST:.c=.o )
OBJ_CAPILLARY   = $(C_CAPILLARY:.c=.o )

LIB_NLS     = $(LIB_MPI) $(LIB_FFT) $(LIB_ADD)

INC_NLS     = $(INC_MPI) $(INC_FFT)
#-------------------------------------------------

# Compilation rule for the executable
# THERE MUST BE A Tab BEFORE THE $(CC) SYMBOL !!!

default:
	@echo "\n   options:   ic, nls, dse, nlstest, nlspost, capillary, hostname\n"

.f.o:
	$(FC) $(FFLAGS) -c $<

nls:
	$(CC)   $(OPT) $(DEF_NLS) $(INC_NLS) -c $(C_NLS)
	$(LINK) $(OPT) $(OBJ_NLS) $(LIB_NLS) -o $(NLS)

dse:
	$(CC)   $(OPT) $(DEF_NLS) $(INC_NLS) -c $(C_DSE)
	$(LINK) $(OPT) $(OBJ_DSE) $(LIB_NLS) -o $(DSE)

ic:
	$(CC)   $(OPT) $(DEF_NLS) $(INC_NLS) -c $(C_IC)
	$(LINK) $(OPT) $(OBJ_IC) $(LIB_NLS) -o $(IC)

nlstest:
	$(CC)   $(OPT) $(DEF_NLS) $(INC_NLS) -c $(C_NLSTEST)
	$(LINK) $(OPT) $(OBJ_NLSTEST) $(LIB_NLS) -o $(NLSTEST)

nlspost:
	$(CC)   $(OPT) $(DEF_NLS) $(INC_NLS) -c $(C_NLSPOST)
	$(LINK) $(OPT) $(OBJ_NLSPOST) $(LIB_NLS) -o $(NLSPOST)

capillary:
	$(CC)   $(OPT) $(DEF_NLS) $(INC_NLS) -c $(C_CAPILLARY)
	$(LINK) $(OPT) $(OBJ_CAPILLARY) $(LIB_NLS) -o $(CAPILLARY)


hostname:
	@echo $(HOSTNAME) 

clean:
	@echo "cleaning ..."
	rm -f *~ *.o

#-------------------------------------------------
