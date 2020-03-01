function twomode_wrap(fbase, gam, ksi)

fnums = 0;
seed  = -1;

gamma   =  [-gam, -gam, 0, 0];    % decay
theta   =  [0, 0, ksi, ksi];      % noise amplitude

n1 = 8;                     % number of timesteps per t=pi, must be multiple of n2
n2 = 8;                     % number of random inputs per t=pi, must be multiple of n3
n3 = 8;                     % number of data saves per t=pi
n4 = 1000;                  % number of pi-periods to run

n=[n1, n2, n3, n4];

%---------------------

lsode_options ("relative tolerance", 1e-10);
lsode_options ("absolute tolerance", 1e-10);
lsode_options ("integration method", "adams");


twomode_core(fbase, seed, gamma, theta, n);

for fnum = 0:fnums

    twomode_core(fbase, fnum, gamma, theta, n);

end

%---------------------

end


