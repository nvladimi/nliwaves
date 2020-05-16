
function twomode_wrap(fbase, G1, G2, F1, F2, dt, isave, nsave, fnums)
%
% "twomode_wrap" is a wrapper script for production runs.
%
% For testing and debugging use "twomode_test" instead.
% Modify for a specific setup and use with "twomode_call" to pass arguments from shell scripts.
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
  
showplot = 1;
seed  = 0;

gamma   =  [-G1, -G1, -G2, -G2];    % decay
force   =  [F1, F1, F2, F2];        % noise amplitude


%---------------------

lsode_options ("relative tolerance", 1e-10);
lsode_options ("absolute tolerance", 1e-10);
lsode_options ("integration method", "adams");


twomode_core(fbase, seed, gamma, force, dt, isave, nsave, showplot);

for fnum = 1:fnums

	     twomode_core(fbase, fnum, gamma, force, dt, isave, nsave, showplot);

end

%---------------------

end


