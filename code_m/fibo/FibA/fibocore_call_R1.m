% This script passes arguments from shell scripts to "fibocore" for production runs.
% input:   fexec=fibocore_call.m;  fout=fbaseout; typ=320;  seed=-1;  
% usage:   octave -qf $fexec  $fout  $typ  $seed  & 

arg_list = argv();

fbase    = arg_list{1};
v        = str2num(arg_list{2});
runtype  = str2num(arg_list{3});
G1       = str2num(arg_list{4});
G2       = str2num(arg_list{5});
dt       = str2num(arg_list{6});
fnum1    = str2num(arg_list{7});

%---------------------------------

isave    = 10;
nsave    = 10000;
fnums    = 5000;
showplot = 0;

m        = 60;
%v        = 1;

P  = 67.65;
%runtype = 320;

%---------------------------------

lsode_options ("relative tolerance", 1e-10);
lsode_options ("absolute tolerance", 1e-10);
lsode_options ("integration method", "adams");

if fnum1 < 0
    seed = fnum1
    fibocore(fbase, seed, P, G1, G2, m, v, runtype, dt, isave, nsave, showplot);
    fnum1 = 1
end

for fnum = fnum1:fnums

      fibocore(fbase, fnum, P, G1, G2, m, v, runtype, dt, isave, nsave, showplot);

end


%--------------------------------
