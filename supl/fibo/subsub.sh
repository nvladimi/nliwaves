#!/bin/bash


#sbatch fibocore_call_r1v38.sub r1_v38L06_dt3 0.6 0.001 -1
#sbatch fibocore_call_r1v58.sub r1_v58R19_dt3 1.9 0.001 -1

sbatch fibocore_call_r1m10.sub r1_m10_R180  0.0 1.80 -1
sbatch fibocore_call_r1m10.sub r1_m10_R190  0.0 1.90 -1

sbatch fibocore_call_r1m50.sub r1_m50_L180  1.80 0.0 -1
sbatch fibocore_call_r1m50.sub r1_m50_L190  1.90 0.0 -1



