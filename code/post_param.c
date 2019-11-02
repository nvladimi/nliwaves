#include "header.h"
#include <time.h>

static int    np, myid;


/*-----------------------------------------------------------*/

void post_param(int *argc, char ***argv, post_ptr post, geom_ptr geom)
{
  char       line[200], param[200], value[200];
  char       filename[80];  
  FILE      *thefile;
  double     pi;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);


  if(argc[0] != 2) {
    if(!myid){printf("\n\t usage: %s run_name\n\n", argv[0][0]);}
    MPI_Finalize();
    exit(1);
  }

  strcpy(post->runname, argv[0][1]);
  strcpy(filename, post->runname);
  strcat(filename, ".post.in");


  /*-- master processor reads the input file and broadcast the parameters --*/

  if (myid==0){

    /* reading input file */

    thefile = fopen (filename, "rt");

    if (thefile == NULL) {
      printf ("\n  Cannot open \"%s\".\n\n", filename);
      MPI_Finalize();
      exit(1);
    }
 
    post->nblur = 0;

    while(fgets(line, 200, thefile) != NULL){

      sscanf(line, "%s  %s", param, value);

      if( strcmp(param,"N0") == 0 )	  geom->N0            = atoi(value);
      if( strcmp(param,"Ngrids") == 0 )	  geom->Ngrids        = atoi(value);
      if( strcmp(param,"Lx") == 0 )	  geom->Lx            = atof(value);
      if( strcmp(param,"Ly") == 0 )	  geom->Ly            = atof(value);
 
      if( strcmp(param,"fn_beg") == 0 )   post->fn_beg       = atoi(value);
      if( strcmp(param,"fn_end") == 0 )   post->fn_end       = atoi(value);
      if( strcmp(param,"fn_inc") == 0 )   post->fn_inc       = atoi(value);
      if( strcmp(param,"bin_out") == 0 )  post->bin_out      = atoi(value);
      if( strcmp(param,"dtime") == 0 )    post->dtime        = atof(value);

      if( strcmp(param,"type") == 0 )     strcpy(post->type, value);
      if( strcmp(param,"rblur") == 0 )    post->rblur[post->nblur++]  = atof(value);
 
      if( strcmp(param,"nbumps") == 0 )   post->nbumps       = atoi(value);
      if( strcmp(param,"rbumps") == 0 )   post->rbumps       = atof(value);
      if( strcmp(param,"hbumps") == 0 )   post->hbumps       = atof(value);

    }

    fclose(thefile);


    /* derived parameters */

    pi =  acos(0)*2;

    if (geom->Lx < 0) geom->Lx = -pi*geom->Lx;
    if (geom->Ly < 0) geom->Ly = -pi*geom->Ly;

    post->time_slices  = (post->fn_end - post->fn_beg + 1) / post->fn_inc;

  }

  /* broadcast the parametrers */
  
  MPI_Bcast( geom, sizeof(geom_str), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast( post, sizeof(post_str), MPI_BYTE, 0, MPI_COMM_WORLD);


  /* checking the geometry */

  if(geom->Lx != geom->Ly) {
    if(!myid){
      printf("\n\t Error: Lx = %f and Ly = %f are not the same. Exiting.\n", 
	     geom->Lx, geom->Ly);
    }
    MPI_Finalize();
    exit(1);
  }


}

/*-----------------------------------------------------------*/
