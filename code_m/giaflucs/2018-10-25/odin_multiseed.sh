#!/bin/bash

tasks=20                        # number of tasks to run simultaneously 
fbase=multiseed                 # base for the file names; parameters in "multiseed.param" 
n1=1                            # first seed
n0=10                           # seeds per task

#-- do not edit below this line ---------------------

i=0
while [ $i -lt $tasks ]
do
     numactl -C +$i octave -qf run_mseed.m $fbase $n1 $n0 &
     let i=i+1
     let n1=n1+n0
done

#wait # IMPORTANT: wait for all to finish or get killed



