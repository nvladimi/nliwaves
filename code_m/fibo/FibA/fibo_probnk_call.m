# This script passes arguments from shell scripts to "fibo_timeavg" for production runs.
# 

arg_list = argv ();

fbase = arg_list{1};
fout  = arg_list{2};

fnum0 = str2num(arg_list{3});

fibo_probnk(fbase, fout, fnum0)


