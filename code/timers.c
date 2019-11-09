#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

#define ON    123
#define OFF  -123
#define NANO  1000000000
#define MAXTIMERS  32

typedef struct
{
  char   name[80];
  int    flag;
  int    ncalls;
  struct timespec dt_tot, dt_min, dt_max, t_on;
} timer_str;


/*---------------------------------------------------------*/

struct timespec diff(struct timespec t1, struct timespec t2);
struct timespec  sum(struct timespec t1, struct timespec t2);
struct timespec  min(struct timespec t1, struct timespec t2);
struct timespec  max(struct timespec t1, struct timespec t2);

extern void timers_resolution_test();
extern void timers_init(char *filebase, double simtime);
extern void timers_out(double simtime);
extern void timers_reset(double simtime);

extern int  timer_set(char *name);
extern void timer_reset(int id);
extern void timer_on(int id);
extern void timer_off(int id);

/*---------------------------------------------------------*/

const struct timespec tzero = {0,0};
const struct timespec tbig  = {99999999,0};

static timer_str     TDB[MAXTIMERS];
static int           ntimers;
static int           myid, np;
static char          filename[80];
static FILE          *thefile;

/*---------------------------------------------------------*/

/**/
void timers_resolution_test()
{
  struct timespec time1, time2, time3, time4;

  printf("\n");
 
  if ( ! clock_getres(CLOCK_REALTIME, &time1))
    printf("REALTIME  resolution   %d:%d\n", (int)time1.tv_sec, (int)time1.tv_nsec);
 
  if (! clock_getres(CLOCK_MONOTONIC, &time2))
    printf("MONOTONIC resolution   %d:%d\n", (int)time2.tv_sec, (int)time2.tv_nsec);

  if (! clock_getres(CLOCK_PROCESS_CPUTIME_ID, &time3))
    printf("PROCESS   resolution   %d:%d\n", (int)time3.tv_sec, (int)time3.tv_nsec);

  if (! clock_getres(CLOCK_THREAD_CPUTIME_ID, &time4))
    printf("THREAD    resolution   %d:%d\n", (int)time4.tv_sec, (int)time4.tv_nsec);

  printf("\n");

}
/**/

/*---------------------------------------------------------*/

void timers_init(char *filebase, double simtime){

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  strcpy(filename, filebase); 
  strcat(filename, ".timers");

  ntimers = 0;

  if (myid == 0) {
    thefile = fopen(filename, "a");
    fprintf(thefile, 
	    "\n# Timers initialized at simulation time t=%12.8e.\n", simtime); 
    fclose(thefile);
  }

}

/*---------------------------------------------------------*/

void timers_reset(double simtime){

  int id;

  for (id=0; id<ntimers; id++) timer_reset(id);

  if (myid == 0) {
    thefile = fopen(filename, "a");
    fprintf(thefile, 
	    "\n# Timers reset at simulation time t=%12.8e.\n", simtime); 
    fclose(thefile);
  }

}

/*---------------------------------------------------------*/

int timer_set(char *name){

  int id;

  id = ntimers++;

  if (id >= MAXTIMERS){
    printf("Max number of timers %d exceeded\n", MAXTIMERS);
    exit(0);
  }

  strcpy(TDB[id].name, name); 

  TDB[id].flag   = OFF;
  TDB[id].dt_tot = tzero;
  TDB[id].dt_min = tbig;
  TDB[id].dt_max = tzero;
  TDB[id].t_on   = tzero;
 
  return(id);

}

/*---------------------------------------------------------*/

void timer_reset(int id){

  TDB[id].flag   = OFF;
  TDB[id].dt_tot = tzero;
  TDB[id].dt_min = tbig;
  TDB[id].dt_max = tzero;
  TDB[id].t_on   = tzero;

}

/*---------------------------------------------------------*/

void timer_on(int id){

  if (TDB[id].flag == ON) {
    printf("An attempt to turn on active timer \"%s\"\n", TDB[id].name);
    exit(0);
  }

  TDB[id].flag = ON;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &TDB[id].t_on);

}

/*---------------------------------------------------------*/

void timer_off(int id){

  struct timespec t_off, dt;

  if (TDB[id].flag == OFF) {
    printf("An attempt to turn off non-active timer \"%s\"\n", TDB[id].name);
    exit(0);
  }

  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t_off);

  dt = diff(TDB[id].t_on, t_off);
  
  TDB[id].dt_tot = sum(TDB[id].dt_tot, dt);
  TDB[id].dt_min = min(TDB[id].dt_min, dt);
  TDB[id].dt_max = max(TDB[id].dt_max, dt);

  TDB[id].ncalls++;
 
  TDB[id].flag   = OFF;

}

/*---------------------------------------------------------*/
 
struct timespec diff(struct timespec t1, struct timespec t2)
{
  struct timespec t;

  if ((t2.tv_nsec - t1.tv_nsec)<0) {
    t.tv_sec  = t2.tv_sec - t1.tv_sec-1;
    t.tv_nsec = NANO + t2.tv_nsec - t1.tv_nsec;
  } else {
    t.tv_sec  = t2.tv_sec  - t1.tv_sec;
    t.tv_nsec = t2.tv_nsec - t1.tv_nsec;
  }

  return t;
}

/*---------------------------------------------------------*/
 
struct timespec sum(struct timespec t1, struct timespec t2)
{
  struct timespec t;
  long int nsec;

  nsec = t2.tv_nsec + t1.tv_nsec;

  if (nsec >= NANO) {
    t.tv_sec  = t2.tv_sec  + t1.tv_sec + 1;
    t.tv_nsec = nsec - NANO;
  } else {
    t.tv_sec  = t2.tv_sec  + t1.tv_sec;
    t.tv_nsec = nsec;
  }

  return t;
}

/*---------------------------------------------------------*/
 
struct timespec min(struct timespec t1, struct timespec t2)
{
  if       (t1.tv_sec < t2.tv_sec) return(t1);
  else if  (t2.tv_sec < t1.tv_sec) return(t2);
  else if  (t1.tv_nsec < t2.tv_nsec) return(t1);
  else return(t2);

}

/*---------------------------------------------------------*/
 
struct timespec max(struct timespec t1, struct timespec t2)
{
  if       (t1.tv_sec > t2.tv_sec) return(t1);
  else if  (t2.tv_sec > t1.tv_sec) return(t2);
  else if  (t1.tv_nsec > t2.tv_nsec) return(t1);
  else return(t2);

}


/*---------------------------------------------------------*/

double tconvert(struct timespec t){

  return( (double)t.tv_sec + (double)t.tv_nsec /NANO );

}
/*---------------------------------------------------------*/

void timers_out(double simtime) {

  int id;

  int    n, n_min, n_max, n_tot; 
  double t, t_tot, t_min, t_max, t_call, t_proc, n_proc;

  /*-- header --*/

  if (myid == 0) {

    thefile = fopen(filename, "a");
    fprintf(thefile, 
	  "\n# Timers values at simulation time t=%12.8e.\n", simtime); 

    fprintf(thefile,"#1.id      2.name       3.t_all_proc 4.t_per_proc   ");
    fprintf(thefile,"5.t_call_min 6.t_call_max 7.t_call_avg   ");
    fprintf(thefile,"8.n_proc_min 9.n_proc_max 10.n_proc_avg\n\n");
 
 }


  /*-- loop over all timers --*/

  for (id=0; id<ntimers; id++){

    /*-- number of calls per process --*/

    n = TDB[id].ncalls;

    MPI_Allreduce(&n, &n_tot, 1,  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (n_tot == 0) continue;

    MPI_Reduce(&n, &n_min, 1,  MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&n, &n_max, 1,  MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);


    n_proc = n_tot/np;


    /*-- times per call --*/ 

    t = tconvert(TDB[id].dt_max);
    MPI_Reduce(&t, &t_max, 1,  MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    t = tconvert(TDB[id].dt_min);
    MPI_Reduce(&t, &t_min, 1,  MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    t = tconvert(TDB[id].dt_tot);
    MPI_Reduce(&t, &t_tot, 1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    t_call = t_tot/n_tot;
    t_proc = t_tot/np;

    /*-- output --*/

    if (myid == 0)
      fprintf(thefile,
	     " %2d  %-16s   %12.6e %12.6e   %12.6e %12.6e %12.6e   %8d %8d %12.2f\n",
	     id, TDB[id].name, t_tot, t_proc, t_min, t_max, t_call,
	     n_min, n_max, n_proc);
  }

  if (myid == 0) fclose(thefile);


}


/*---------------------------------------------------------*/
