#include "header.h"
#include <time.h>

static int    np, myid;
       char   msg[80];

void print_header(geom_ptr geom, ctrl_ptr ctrl, 
                  phys_ptr phys, diag_ptr diag);


/*-----------------------------------------------------------*/

void parameters(int *argc, char ***argv, 
		geom_ptr geom, ctrl_ptr ctrl, 
                phys_ptr phys, diag_ptr diag)
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

  /*-- zero out structures --*/

  memset(geom, 0, sizeof(*geom));
  memset(ctrl, 0, sizeof(*ctrl));
  memset(phys, 0, sizeof(*phys));
  memset(diag, 0, sizeof(*diag));


  /*-- get current run name and the name of parameter file --*/

  strcpy(ctrl->runname, argv[0][1]);
  strcpy(filename, ctrl->runname);
  strcat(filename, ".in");

 
  /*-- master processor reads the input file and broadcasts the parameters --*/

  if (myid==0){

    /*-- reading input file --*/

    thefile = fopen (filename, "rt");

    if (thefile == NULL) {
      printf ("\n  Cannot open \"%s\".\n\n", filename);
      MPI_Finalize();
      exit(1);
    }
 
    strcpy(phys->f_type, "none");

    while(fgets(line, 200, thefile) != NULL){

      sscanf(line, "%s  %s", param, value);

      if( strcmp(param,"Ngrids") == 0 )	  geom->Ngrids        = atoi(value);
      if( strcmp(param,"N0") == 0 )	  geom->N0            = atoi(value);
      if( strcmp(param,"Lx") == 0 )	  geom->Lx            = atof(value);
      if( strcmp(param,"Ly") == 0 )	  geom->Ly            = atof(value);
      if( strcmp(param,"dealiasZ") == 0 ) geom->dealiasZ      = atof(value);
      if( strcmp(param,"dealiasF") == 0 ) geom->dealiasF      = atoi(value);
      if( strcmp(param,"regridZup") == 0 ) geom->regridZup    = atof(value);
      if( strcmp(param,"regridZdn") == 0 ) geom->regridZdn    = atof(value);
      if( strcmp(param,"regridTh") == 0 )  geom->regridTh     = atof(value);
      if( strcmp(param,"sieve")    == 0 )  geom->sieve        = atoi(value);
      if( strcmp(param,"klow")    == 0 )   geom->klow         = atoi(value);

      if( strcmp(param,"coefA") == 0 )    phys->coefA         = atof(value);
      if( strcmp(param,"coefB") == 0 )    phys->coefB         = atof(value);
      if( strcmp(param,"coefC") == 0 )    phys->coefC         = atof(value);
      if( strcmp(param,"coefE") == 0 )    phys->coefE         = atof(value);
      if( strcmp(param,"coefNu") == 0 )   phys->coefNu        = atof(value);
      if( strcmp(param,"coefRho") == 0 )  phys->coefRho       = atof(value);
      if( strcmp(param,"expoS") == 0 )    phys->expoS         = atof(value);
      if( strcmp(param,"focus") == 0 )    phys->focus         = atoi(value);
      if( strcmp(param,"mask") == 0 )     phys->mask          = atoi(value);

      if( strcmp(param,"f_type") == 0 )   strcpy(phys->f_type, value);

      if( strcmp(param,"f_kmin") == 0 )   phys->f_kmin        = atof(value);
      if( strcmp(param,"f_kmax") == 0 )   phys->f_kmax        = atof(value);
      if( strcmp(param,"f_kdamp") == 0 )  phys->f_kdamp       = atof(value);
      if( strcmp(param,"f_kfric") == 0 )  phys->f_kfric       = atof(value);
      if( strcmp(param,"f_alpha") == 0 )  phys->f_alpha       = atof(value);
      if( strcmp(param,"f_beta") == 0 )   phys->f_beta        = atof(value);
      if( strcmp(param,"f_radius") == 0 ) phys->f_radius      = atof(value);


      if( strcmp(param,"tmax")    == 0 )  ctrl->tmax          = atof(value);
      if( strcmp(param,"dnInfo")  == 0 )  ctrl->dnInfo        = atoi(value);
      if( strcmp(param,"dnDiag")  == 0 )  ctrl->dnDiag        = atoi(value);
      if( strcmp(param,"dnData")  == 0 )  ctrl->dnData        = atoi(value);
      if( strcmp(param,"njig")    == 0 )  ctrl->njig          = atoi(value);
      if( strcmp(param,"noff")    == 0 )  ctrl->noff          = atoi(value);
      if( strcmp(param,"dt0")     == 0 )  ctrl->dt0           = atof(value);
      if( strcmp(param,"CFL")     == 0 )  ctrl->CFL           = atof(value);

      if( strcmp(param,"pdfMax")   == 0 ) diag->pdfMax        = atof(value);
      if( strcmp(param,"pdfNbins") == 0 ) diag->pdfNbins      = atoi(value);
      if( strcmp(param,"clpsNmax") == 0 ) diag->clpsNmax      = atoi(value);
      if( strcmp(param,"clpsCmax") == 0 ) diag->clpsCmax      = atoi(value);
      if( strcmp(param,"clpsPsi0") == 0 ) diag->clpsPsi0      = atof(value);
      if( strcmp(param,"clpsRate") == 0 ) diag->clpsRate      = atof(value);
      if( strcmp(param,"vortexTh") == 0 ) diag->vortexTh      = atof(value);

    }

    fclose(thefile);


    /*-- set restart time and file number from tag file --*/

    sprintf(filename,"%s.tag", ctrl->runname);
    thefile = fopen (filename, "rt");
    if (thefile != NULL) {
      fscanf(thefile, "%d  %d  %d\n", 
	     &(ctrl->restart), &(ctrl->n0), &(geom->grid));
      fclose(thefile); 
    } else {
      sprintf(msg, "#err:  tag file not found. Exiting.\n");
      io_message(msg,0);
      exit(-1);
    }

    if ( geom->grid >=  geom->Ngrids ) {
      sprintf(msg, "#err:  insufficent number of grids. Exiting.\n");
      io_message(msg,0);
      exit(-1);
    }

    /*-- derived parameters --*/

    pi =  acos(0)*2;

    if (geom->Lx < 0) geom->Lx = -pi*geom->Lx;
    if (geom->Ly < 0) geom->Ly = -pi*geom->Ly;

    print_header(geom, ctrl, phys, diag);

  }

  /*-- broadcast the parametrers --*/
  
  MPI_Bcast( geom, sizeof(geom_str), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast( ctrl, sizeof(ctrl_str), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast( phys, sizeof(phys_str), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast( diag, sizeof(diag_str), MPI_BYTE, 0, MPI_COMM_WORLD);


  /*-- checking the geometry --*/

  if(geom->Lx != geom->Ly) {
    sprintf(msg,
	    "\n\t Error: Lx = %f and Ly = %f are not the same. Exiting.\n", 
	    geom->Lx, geom->Ly);
    io_message(msg,0);
    exit(-1);
  }


}



/*-----------------------------------------------------------*/

void print_header(geom_ptr geom, ctrl_ptr ctrl, 
		  phys_ptr phys, diag_ptr diag)
{

  char       filename[80];  
  FILE      *thefile;
  time_t     tmr;
  char*      date;
 
  if (myid != 0) return;

  /*-- get current time --*/

  tmr = time(NULL);
  date = ctime(&tmr); 


  /*-- write text information --*/

  strcpy (filename, ctrl->runname);
  strcat (filename, ".out");  thefile = fopen (filename, "a");

  fprintf(thefile, "\n");
  fprintf(thefile, "# Started:   %s", date);  
  fprintf(thefile, "# Processes: %d\n", np);  
  fprintf(thefile, "# Run name:  %s\n", ctrl->runname );

  fprintf(thefile, "#\n# Domain geometry: \n");
  fprintf(thefile, "#     Ngrids= %d\n",  geom->Ngrids);
  fprintf(thefile, "#     N0    = %d\n",  geom->N0);
  fprintf(thefile, "#     Lx    = %f\n",  geom->Lx);
  fprintf(thefile, "#     Ly    = %f\n",  geom->Ly);
  fprintf(thefile, "#     dealiasZ  = %f\n",  geom->dealiasZ);
  fprintf(thefile, "#     dealiasF  = %d\n",  geom->dealiasF);
  fprintf(thefile, "#     regridZup = %f\n",  geom->regridZup);
  fprintf(thefile, "#     regridZdn = %f\n",  geom->regridZdn);
  fprintf(thefile, "#     regridTh  = %e\n",  geom->regridTh);
  fprintf(thefile, "#     sieve    = %d\n",   geom->sieve);
  fprintf(thefile, "#     klow     = %d\n",   geom->klow);

  fprintf(thefile, "#\n# Physics parametes: \n");
  fprintf(thefile, "#     coefA   = %11.4e \n", phys->coefA);
  fprintf(thefile, "#     coefB   = %11.4e \n", phys->coefB);
  fprintf(thefile, "#     coefC   = %11.4e \n", phys->coefC);
  fprintf(thefile, "#     coefE   = %11.4e \n", phys->coefE);
  fprintf(thefile, "#     coefNu  = %11.4e \n", phys->coefNu);
  fprintf(thefile, "#     coefRho = %11.4e \n", phys->coefRho);
  fprintf(thefile, "#     expoS   = %11.4e \n", phys->expoS);
  fprintf(thefile, "#     focus   = %3d\n",     phys->focus);
  fprintf(thefile, "#     mask    = %3d\n",     phys->mask);

  fprintf(thefile, "#     f_type  = %s\n",    phys->f_type);
  fprintf(thefile, "#     f_kmin  = %f\n",    phys->f_kmin);
  fprintf(thefile, "#     f_kmax  = %f\n",    phys->f_kmax);
  fprintf(thefile, "#     f_kdamp = %f\n",    phys->f_kdamp);
  fprintf(thefile, "#     f_kfric = %f\n",    phys->f_kfric);
  fprintf(thefile, "#     f_alpha = %f\n",    phys->f_alpha);
  fprintf(thefile, "#     f_beta  = %f\n",    phys->f_beta);
  fprintf(thefile, "#     f_radius= %f\n",    phys->f_radius);

  fprintf(thefile, "#\n# Input and output:\n");
  fprintf(thefile, "#     tmax     = %f\n", ctrl->tmax);
  fprintf(thefile, "#     dnInfo   = %d\n", ctrl->dnInfo);
  fprintf(thefile, "#     dnDiag   = %d\n", ctrl->dnDiag);
  fprintf(thefile, "#     dnData   = %d\n", ctrl->dnData);
  fprintf(thefile, "#     njig     = %d\n", ctrl->njig);
  fprintf(thefile, "#     noff     = %d\n", ctrl->noff);
  fprintf(thefile, "#     dt0      = %f\n", ctrl->dt0);
  fprintf(thefile, "#     CFL      = %f\n", ctrl->CFL);

  fprintf(thefile, "#\n# Diagnostics:\n");
  fprintf(thefile, "#     pdfMax   = %f\n", diag->pdfMax);
  fprintf(thefile, "#     pdfNbins = %d\n", diag->pdfNbins);
  fprintf(thefile, "#     clpsNmax = %d\n", diag->clpsNmax);
  fprintf(thefile, "#     clpsCmax = %d\n", diag->clpsCmax);
  fprintf(thefile, "#     clpsPsi0 = %f\n", diag->clpsPsi0);
  fprintf(thefile, "#     clpsRate = %f\n", diag->clpsRate);
  fprintf(thefile, "#     vortexTh = %f\n", diag->vortexTh);

  fprintf(thefile, "#\n# Restart:\n");
  fprintf(thefile, "#     restart  = %d\n", ctrl->restart);
  fprintf(thefile, "#     n0       = %d\n", ctrl->n0);
  fprintf(thefile, "#     grid     = %d\n", geom->grid);

  fprintf(thefile, "\n");


  fclose(thefile);

}

/*-----------------------------------------------------------*/

