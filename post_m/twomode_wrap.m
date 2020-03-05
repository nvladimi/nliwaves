function twomode_wrap(fbase, gam, ksi, dt, isave, nsave)

fnums = 0;
seed  = -1;

gamma   =  [-gam, -gam, 0, 0];    % decay
theta   =  [0, 0, ksi, ksi];      % noise amplitude


%---------------------

lsode_options ("relative tolerance", 1e-10);
lsode_options ("absolute tolerance", 1e-10);
lsode_options ("integration method", "adams");


twomode_core(fbase, seed, gamma, theta, dt, isave, nsave);

for fnum = 0:fnums

	     twomode_core(fbase, fnum, gamma, theta, dt, isave, nsave);

end

%---------------------

end


