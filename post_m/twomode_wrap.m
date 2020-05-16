function twomode_wrap(fbase, gam, ksi, dt, isave, nsave, fnums)

seed  = -1;

gamma   =  [0, 0, -gam, -gam];    % decay
theta   =  [ksi, ksi, 0, 0];      % noise amplitude


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


