% This script passes arguments from shell scripts to "fibocore" for production runs.

% for use in multiple runs add input arguments, for example:
%
% arg_list = argv();
% fbase     = arg_list{1};
% alpha     = str2num(arg_list{2});

%---------------------------------


fbase="dblB/dblB";

runtype="dblB";

alpha=0.375;
dt=0.001;
G1=0.1;
G2=0;
fnum1=-1;


isave    = 10;
nsave    = 10;
fnums    = 5;
showplot = 0;

m        = 60;
mp       = 50;

PI       = 10.0;

%---------------------------------

lsode_options ("relative tolerance", 1e-10);
lsode_options ("absolute tolerance", 1e-10);
lsode_options ("integration method", "adams");

if fnum1 < 0
    seed = fnum1;
    fibocore(runtype, fbase, seed, alpha, m, mp, PI, G1, G2, dt, isave, nsave, showplot);
    fnum1 = 1;
end

for fnum = fnum1:fnums
    fibocore(runtype, fbase, fnum, alpha, m, mp, PI, G1, G2, dt, isave, nsave, showplot);
end


%--------------------------------
