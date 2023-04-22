#!/bin/bash

#--  Before running, create a parameter file, "set001.param", etc, for each set. 
#--  Submitting the job, use number of tasks equal to the number of sets.
 
tasks=3

#-- do not edit below this line -------------------

i=0
n=1
while [ $i -lt $tasks ]
do
     fbase=$(printf set%03d $n)
     numactl -C +$i octave -qf run_mset.m $fbase &
     let i=i+1
     let n=n+1
done

#wait # IMPORTANT: wait for all to finish or get killed

