/*---

  gcc test_timers.c timers.c -I/usr/lib/mpich/include  -lrt -lmpi 

 ---*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern void timers_resolution_test();
extern void timers_init(char *filebase, double simtime);
extern void timers_out(double simtime);
extern void timers_reset(double simtime);

extern int  timer_set(char *name);
extern void timer_reset(int id);
extern void timer_on(int id);
extern void timer_off(int id);


/*------------------------------------------------------------*/
 
int main(int argc, char **argv)
{

  int temp, i,k;
  int tmr1, tmr2;

  MPI_Init(&argc, &argv);

  timers_resolution_test();

  timers_init("timers_test", 0);


  tmr1= timer_set("outer_loop");
  tmr2= timer_set("inner_loop");

  timer_on(tmr1);

  for (k=0; k<4; k++) {

    timer_on(tmr2);
    for (i=0; i<100000000; i++) temp+=temp;
    timer_off(tmr2);
  }

  timer_off(tmr1);

  timers_out(1);

  MPI_Finalize();

  return 0;
}
 
