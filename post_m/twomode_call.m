
arg_list = argv ();

fbase = arg_list{1};
gam   = str2num(arg_list{2});
ksi   = str2num(arg_list{3});
dt    = str2num(arg_list{4});
isave = str2num(arg_list{5});
nsave = str2num(arg_list{6});
fnums = str2num(arg_list{7});

twomode_wrap(fbase, gam, ksi, dt, isave, nsave, fnums)


