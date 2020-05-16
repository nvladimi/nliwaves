
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
% gamma        strength of decay (<0) or multiplicative forcing (>0), gamma[re1, im1, re2, im2]
% force        noise amplitude, force[re1, im1, re2, im2]
% dt           timestep
% isave        save data every "isave" timestep
% nsave        generate "nsave" number of saves
% showplot     supress debudding plots if showplot=0
%

fbase     = "tm0";
seed     = 0;
gamma    = [0, 0, 0, 0];
force    = [0, 0, 0, 0]; 
dt       = 0.4;
isave    = 1;
nsave    = 100;
showplot = 1;

%---------------------

lsode_options ("relative tolerance", 1e-10);
lsode_options ("absolute tolerance", 1e-10);
lsode_options ("integration method", "adams");


twomode_core(fbase, seed, gamma, force, dt, isave, nsave, showplot);

%---------------------

end


