#include <stdio.h>
#include <fftw_mpi.h>
#include <mpi.h>

#define space_loop for(x=0;x<lnx;++x)for(y=0;y<ny;++y)
#define space_loop2 for(x=0;x<1;++x)for(y=0;y<2;++y)

/* 
 *  Athelas:
 *    gcc nersc_example.c  -I/usr/lib/mpich/include -lmpi -lm -lfftw_mpi
 *
 *  Coyote:
 *      module load fftw/2.1.5  openmpi-gcc/1.2.4
 *      gcc nersc_example.c -I$FFTW_INCLUDE -I$MPIHOME/include 
 *        -L$MPIHOME/lib64 -lmpi -lm -L$FFTW_HOME/lib -lfftw_mpi -lfftw
 */

/*
 * 	nx*ny*nz FFT computed in parallel (MPI FFTW)
 *	
 *	slab decomposed (along x) among numtasks MPI tasks
 *		
 *	lnx, lny, lny		 local nx,ny,nz
 *	lnyt			 local ny after transpose
 *	lxs			 local x start
 *	lyst			 local y start after transpose
 *	lsize			 total local size
 *		
 *	Frank Hale Jul 2008 (fvhale@nersc.gov)
 *	David Skinner	Aug 2001 (dskinner@nersc.gov)
 */

int main(int argc, char **argv)
{
     int 		nn, nx,ny; 		        /* DIMS */
     int 		myid; 			        /* MPI */
     int 		lnx, lxs, lnyt, lyst, lsize;	/* LOCAL DIMS */
     int 		x, y;   			/* LOOP VARS */
     fftwnd_mpi_plan 	plan, iplan;			/* PLANS */
     fftw_complex 	*data, *work;			/* VECTORS */

     MPI_Init(&argc,&argv);

     if(argc != 2) {
      if(!myid){printf("usage: %s dim_size\n",argv[0]);}
      if(!myid){printf("\tdim_size = number of points in each dimension\n");}
      MPI_Finalize();
      exit(1);
     }
     else { nn = atoi(argv[1]); }

     MPI_Comm_rank(MPI_COMM_WORLD, &myid);

     nx = ny = nn;
     printf("PE: %d, nx=%d, ny=%d\n", myid, nx, ny);

     plan = fftw2d_mpi_create_plan(MPI_COMM_WORLD, nx, ny,
	 FFTW_FORWARD, FFTW_ESTIMATE);

     iplan = fftw2d_mpi_create_plan(MPI_COMM_WORLD, nx, ny,
	 FFTW_BACKWARD, FFTW_ESTIMATE);

/*   slab decomposition */

     fftwnd_mpi_local_sizes(plan, &lnx, &lxs, &lnyt, &lyst, &lsize);

     printf("data_decomp:: %d xslab %d->%d (%d) %d (%.3e MB of %.3e MB total)\n"
	,myid, lxs, lxs+lnx,lnx, lyst,
	lsize*sizeof(fftw_complex)*1.0e-6,
	nn*nn*sizeof(fftw_complex)*1.0e-6);

/*  malloc data and work vectors */

     data = (fftw_complex*) malloc(sizeof(fftw_complex) * lsize);
     if(!data) {printf("malloc failed for data (task = %d)\n",myid);}

     work = (fftw_complex*) malloc(sizeof(fftw_complex) * lsize);
     if(!data) {printf("malloc failed for work (task = %d)\n",myid);}

/* initialize data */

     space_loop{
        data[x*ny + y][0] = 1.0;
        data[x*ny + y][1] = 0.0;
     }

     space_loop2{
        printf("PE %d input %d %d %.3e %.3e\n", myid, x+lxs, y,
	 data[x*ny + y][0],
	 data[x*ny + y][1]);
     }

     MPI_Barrier(MPI_COMM_WORLD);


     fftwnd_mpi(plan, 1, data, work, FFTW_TRANSPOSED_ORDER);

/* after forward fft */

     space_loop2{
        printf("PE %d forward %d %d %.3e %.3e\n", myid, x+lxs, y,
	 data[x*ny + y][0],
	 data[x*ny + y][1]);
     }

     MPI_Barrier(MPI_COMM_WORLD);
     fftwnd_mpi(iplan, 1, data, work, FFTW_TRANSPOSED_ORDER);

/* after reverse fft */

     space_loop2{
        printf("PE %d reverse %d %d %.3e %.3e\n", myid, x+lxs, y,
	 data[x*ny + y][0],
	 data[x*ny + y][1]);
     }

     MPI_Barrier(MPI_COMM_WORLD);
     printf("\nWrap up PE %d\n",myid);

     free(data);
     free(work);
     fftwnd_mpi_destroy_plan(plan);
     fftwnd_mpi_destroy_plan(iplan);

     MPI_Finalize();
     exit(0); 
}
