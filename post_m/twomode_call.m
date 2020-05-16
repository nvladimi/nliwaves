# This script passes arguments from shell scripts to "twomode_wrap" for production runs.
# 
# debug:  octave -qf twomode_call.m "tm4" 0.1 0  0  0  0.4  1 100  0  
# run:    numactl -C +0  octave -qf twomode_call.m "xi2b/xi2b" 2.5e-2 0 0 1.25e-2  0.1  1 10000 200  &

arg_list = argv ();

fbase = arg_list{1};
G1    = str2num(arg_list{2});
G2    = str2num(arg_list{3});
F1    = str2num(arg_list{4});
F2    = str2num(arg_list{5});
dt    = str2num(arg_list{6});
isave = str2num(arg_list{7});
nsave = str2num(arg_list{8});
fnums = str2num(arg_list{9});

twomode_wrap(fbase, G1, G2, F1, F2, dt, isave, nsave, fnums)


