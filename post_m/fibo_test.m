
function fibo_test
%
% "twomode_test" is a wrapper script for "twomode_core" for test runs.
% For production runs use "twomode_wrap" instead.
%
% Parameters to "twomode_core":
%
% fbase        sring base for output files
% seed/fnum    when positive, second argument to "twomode_core" restore IC from (fnum-1) files
%              otherwise use as seed to create (fnum=0) file
% G            level of damping (application depends on runtype)
% P            level of forcing (application depends on runtype)
% M            number of modes
% runtype      integer to distinguish different forcing and damping and IC
% dt           timestep
% isave        save data every "isave" timestep
% nsave        generate "nsave" number of saves
% showplot     supress debudding plots if showplot=0
%
%  first digit in "runtype" determines initial condition
%  second digit determines forcing
%
%  ic_type 1:        a0 =  [1, 1, 0, ..., 0]
%  ic_type 2:        a0 =  [1, 1, -1/sqrt(2), 0, ..., 0]
%  ic_type 3:        |a0|^2 = 1/Fi,   random phases
%  ic_type 4:        |a0|^2 = random, random phases 
%
%  force_type 0:     force = 0, gamma = 0;
%  force_type 1:     force = [P, P, 0, ..., 0],         gamma = [0, ..., 0, G, G];        direct cascade
%  force_type 2:     force = [P, P, 0, ..., 0],         gamma = [0, ..., 0, G*Fi, G*Fi];  direct cascade balanced  
%  force_type 3:     force = [0, ..., 0, P, P],         gamma = [G, G, 0, ..., 0];        inverse cascade
%  force_type 4:     force = [0, ..., 0, P/Fi, P/Fi],   gamma = [G, G, 0, ..., 0];        inverse cascage balanced
%  force_type 5:     force = [P, P, P, ..., P],         gamma = P*Fi/2;                   equilibrium by force
%  force_type 6:     force = 2*G*Fi,                    gamma = [G, G, G, ..., G];        equilibrium by decay
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
runtype  = 20;

%---------------------

lsode_options ("relative tolerance", 1e-10);
lsode_options ("absolute tolerance", 1e-10);
lsode_options ("integration method", "adams");


fibo_core(fbase, seed, G, P, M, runtype, dt, isave, nsave, showplot);

%---------------------

end


