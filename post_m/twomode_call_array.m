# This script passes arguments from shell scripts to "twomode_array_n1n2theta" for production runs.
# 

arg_list = argv ();

fbase = arg_list{1};
fnum1 = str2num(arg_list{2});
fnum2 = str2num(arg_list{3});
fbaseout = arg_list{4};

fnums = [fnum1:fnum2];


twomode_array_n1n2theta(fbase, fnums, fbaseout)

