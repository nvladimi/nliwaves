
function twomode_test
%
% "twomode_test" is a wrapper script for test runs.
% For production runs use "twomode_wrap" instead.
%
% Parameters to "twomode_core":
%
% fbase        sring base for output files
% seed/fnum    when positive, second argument to "twomode_core" restore IC from (fnum-1) files
%              otherwise use as seed to create (fnum=0) file
% Gamma        strength of decay (<0) or multiplicative forcing (>0), Gamma[re1, im1, re2, im2]
% Rflux        flux of random force, Rflux[re1, im1, re2, im2]
% dt           timestep
% isave        save data every "isave" timestep
% nsave        generate "nsave" number of saves
% showplot     supress debudding plots if showplot=0
%

fbase     = "tm0";
seed     = 0;
dt       = 0.1;
isave    = 1;
nsave    = 1000;
showplot = 1;


%-- Direct -- 
Gamma    = [0, 0, -1, -1]*0.01;
Rflux    = [1, 1, 0, 0] * 4e-3; 

%-- Inverse -- 
Gamma    = [-1, -1, 0, 0]*0.01;
Rflux    = [0, 0, 1, 1] * 1e-3; 

%-- Equilibrium
Gamma    = [-1, -1, -2, -2]*0.01;
Rflux    = [1, 1, 1, 1] * 1e-3; 


%---------------------

lsode_options ("relative tolerance", 1e-10);
lsode_options ("absolute tolerance", 1e-10);
lsode_options ("integration method", "adams");


twomode_core(fbase, seed, Gamma, Rflux, dt, isave, nsave, showplot);

%---------------------

end


