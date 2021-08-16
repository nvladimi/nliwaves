# This script passes arguments from shell scripts to "fibo_timeavg" for production runs.
# 

arg_list = argv ();

fbase = arg_list{1};
fout  = arg_list{2};

fnum1 = str2num(arg_list{3});
fnum2 = str2num(arg_list{4});

fnums = [fnum1:fnum2];


fibo_moments(fbase, fout, fnums)


