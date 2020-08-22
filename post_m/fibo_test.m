
function fibo_test
%
% "fibo_test" is a wrapper script for "fibocore" for test runs.
% For production runs use "fibo_wrap" instead.
%
% Parameters to "twomode_core":
%
% fbase        sring base for output files
% seed/fnum    when positive, second argument to "twomode_core" restore IC from (fnum-1) files
%              otherwise use as seed to create (fnum=0) file
% P            level of forcing of a specified middle mode 
% G1           level of damping at low modes
% G2           level of damping at high modes
% M            number of modes
% runtype      integer to distinguish different IC
% dt           timestep
% isave        save data every "isave" timestep
% nsave        generate "nsave" number of saves
% showplot     supress debudding plots if showplot=0
%
%  first digit in "runtype" determines initial condition
%  second two digits determine forcing mode
%
%  ic_type 1:        a0 =  [1, 1, 0, ..., 0]
%  ic_type 2:        a0 =  [1, 1, -1/sqrt(2), 0, ..., 0]
%  ic_type 3:        |a0|^2 = 1/Fi,   random phases
%  ic_type 4:        |a0|^2 = random [0:1] , random phases 
%  ic_type 5:        a0 = 0
%

fbase    = "f0";
seed     = 0;
dt       = 0.1;
isave    = 1;
nsave    = 200;
showplot = 1;

G        = 0;
P        = 0;
M        = 3;
runtype  = 200;

%---------------------

lsode_options ("relative tolerance", 1e-10);
lsode_options ("absolute tolerance", 1e-10);
lsode_options ("integration method", "adams");


fibo_core(fbase, seed, P, G, G, M, runtype, dt, isave, nsave, showplot);

%---------------------

end


