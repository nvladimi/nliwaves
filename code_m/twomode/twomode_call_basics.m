# This script passes arguments from shell scripts to "twomode_basics" for production runs.
# 

arg_list = argv ();

fbase = arg_list{1};
fnum2 = str2num(arg_list{2});

fnum1 = 1; 
fnums = [fnum1:fnum2];
ev    = 0;

twomode_basics(fbase, fnums, ev)


