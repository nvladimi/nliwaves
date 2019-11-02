function h_from_ode

%-- parameters --

  u0 = [0.38,  -0.67,6.6];     % initial [L, L_t, beta] at t=0.1
  u0 = [0.231, -2.09,0.461];   % initial [L, L_t, beta] at t=0.2
  u0 = [0.12785,-3.343,0.2027];  % initial [L, L_t, beta] at t=0.24

  fname1='test1.dat';
  fname2='test2.dat';

  dt   = 1.e-5;
  tmax = 0.1;

%-- matlab -- different from octave interface to ODE solver
%
%  options = odeset('MaxStep',dt);
%
%  [t,u1]=ode45( @udot_eq21, [dt,tmax], u0, options);
%  [t,u2]=ode45( @udot_eq22, [dt,tmax], u0, options);
%  t=t';
  

%-- octave -- different from matlab interface to ODE solver

   t = dt:dt:tmax;

   %u1 = lsode('udot_eq21', u0, t);
   u2 = lsode('udot_eq22', u0, t);

%-- the rest is similar

   %printfile(fname1, u0, u1, t);
   printfile('eq22guess.dat', u0, u2, t);

end

%---------------------

function du = udot_eq21(u,r)   % octave use udot(u,r);  matlab use udot(r,u)
%
% model given by Eq.21 in the proposal  

  Nc =  11.7008822458640;
  P  =  13.8944704347304;
  M  =  3.47361760868355;
  P  =  18.65;

  epsilon = 0.01;
  a       = 1;

  c1 = 1;
  c2 = 2*epsilon*(a+2) * Nc/M;
  c3 = a/(a+2) * (P-Nc)/Nc;
  c4 = epsilon*( (P-Nc)/M - 2*a );

  %if u(3) > 0 
  %    nu = 45.1*exp( -pi * u(3)^(-0.5) );
  %else
  %    nu = 0;
  %end


  du = [0; 0; 0];
  
  du(1) =  u(2);
  du(2) = -c1* u(3) / u(1)^3;
  du(3) = -c2 * (1 + 0.5*c3*u(3)) / u(1)^2 + c4*u(2)^2; %-nu/u(1)^2 ; 

end

%---------------------

function du = udot_eq22(u,r)   % octave use udot(u,r);  matlab use udot(r,u)
%
% model given by Eq.22 in the proposal  

  epsilon = 0.01;
  V    =  2.20586;
  dBdV = -1.9216523;
  

  c1 = 1 - 8*V/dBdV;
  c2 = -epsilon*V*dBdV;
  c3 = 1;
  c4 = 0;

  c1 = 0.6;
  c2 = 0.2;
  c3 = -0.6;

  du = [0; 0; 0];
  
  du(1) =  u(2);
  du(2) = -c1* u(3) / u(1)^3;
  du(3) = -(c2 + c3*u(3)) / u(1)^2 + c4*u(2)^2;

end

%---------------------

function printfile(fname, u0, u, t)

   h = Vbeta(u(:,3))./u(:,1);

   fid  = fopen(fname, 'wt');
   fprintf(fid, '#_1.t   2.L   3.L_t   4.beta  5.h \n\n');
   fprintf(fid, '%8.4f  %12.4e  %12.4e  %12.4e  %12.4e\n', ...
                0, u0(1), u0(2), u0(3), 0);

   for i = 2:length(t)
     fprintf(fid, '%8.4e  %12.4e  %12.4e  %12.4e  %12.4e\n',...
                   t(i), u(i,1), u(i,2), u(i,3), h(i));
   end

   fclose(fid);

end

%---------------------

function b = betaV(v)

  b6=      98.9922359479098;
  b5=   -1403.7339662171987;
  b4=    8286.5470292400041;
  b3=  -26067.2261230988006;
  b2=   46084.0997339975365;
  b1=  -43410.0271983773928;
  b0=   17021.0686019477034;

  b = b0 + b1*v + b2*v.^2 + b3*v.^3 + b4*v.^4 + b5*v.^5 + b6*v.^6;

end


%---------------------

function v = Vbeta(b)

  a8=  -1.052411266450615;
  a7=  -3.188848891246191;
  a6=  -2.514007386433376;
  a5=   1.005165659608414;
  a4=   2.015195330403522;
  a3=   0.442514242344889;
  a2=  -0.423409879090351;
  a1=  -0.549418771419909;
  a0=   2.203522026963209;

  v = a0 + a1*b + a2*b.^2 + a3*b.^3 + a4*b.^4 + ...
        a5*b.^5 + a6*b.^6 + a7*b.^7 + a8*b.^8;

end

%---------------------
