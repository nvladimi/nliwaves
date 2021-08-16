% This script passes arguments from shell scripts to "fibocore" for production runs.

% for use in multiple runs add input arguments, for example:
%
% arg_list = argv();
% fbase     = arg_list{1};
% alpha     = str2num(arg_list{2});

%---------------------------------


fbase="fibB/fibB"; %frun=g1_v38L06

runtype="fibB";

alpha=0.375;
dt=0.001;
G1=0.6;
G2=0;
fnum1=-1;


isave    = 10;     %10
nsave    = 10;     %10000
fnums    = 5;      %5000
showplot = 0;

m        = 60;
mp       = 50;

PI       = 67.65;

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
