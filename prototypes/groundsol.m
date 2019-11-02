function groundsol

  beta = 0.00;  u01 = 2.206200;  u02 = 2.206201;  fname='beta0.00.dat';

  u01 = 2.20620086465074607478363406;
  NP = 11.700896524559653878655;

  %beta = 0.01;  u01 = 2.200965;  u02 = 2.200966;  fname='beta0.01.dat';
  %beta = 0.02;  u01 = 2.195610;  u02 = 2.195611;  fname='beta0.02.dat';
  %beta = 0.03;  u01 = 2.190122;  u02 = 2.190123;  fname='beta0.03.dat';
  %beta = 0.10;  u01 = 2.14663;   u02 = 2.14664;   fname='beta+0.10.dat';
  %beta = 0.20;  u01 = 2.064;     u02 = 2.065;     fname='beta+0.20.dat';
  %beta = 0.30;  u01 = 1.97;      u02 = 1.98;      fname='beta+0.30.dat';

  %beta = -0.01; u01 = 2.211325;  u02 = 2.211326;  fname='beta-0.01.dat';
  %beta = -0.02; u01 = 2.216347;  u02 = 2.216348;  fname='beta-0.02.dat';
  %beta = -0.03; u01 = 2.221272;  u02 = 2.221273;  fname='beta-0.03.dat';
  %beta = -0.04; u01 = 2.226107;  u02 = 2.226108;  fname='beta-0.04.dat';
  %beta = -0.05; u01 = 2.230856;  u02 = 2.230857;  fname='beta-0.05.dat';
  %beta = -0.10; u01 = 2.253462870;  u02 = 2.253462871;  fname='beta-0.10.dat';
  %beta = -0.20; u01 = 2.294153566;  u02 = 2.294153567;  fname='beta-0.20.dat';
  %beta = -0.30; u01 = 2.330380325;  u02 = 2.330380326;  fname='beta-0.30.dat';

  fname='test.dat';

  dr   = 1.e-3;
  rmax = 20;

%-- matlab -- different from octave interface to ODE solver
%
%  options = odeset('MaxStep',dr);
%
%  [r,u1]=ode45( @udot2d, [dr,rmax], ini(u01,dr,beta), options);
%  %[r,u2]=ode45( @udot2d, [dr,rmax], ini(u02,dr,beta), options);
%  r=r';
  
%-- octave

   r = dr:dr:rmax;

   u1 = lsode('udot2d', ini(u01,dr,beta), r);
   u2 = lsode('udot2d', ini(u02,dr,beta), r);

%-- the rest is similar

   plot(r, u1(:,1));
   %plot(r, u1(:,1), r, u2(:,1), r, 0*r, r, 0.5*sqrt(-beta)*r );
   %semilogy(r, abs(u1(:,1)), r, abs(u2(:,1))); %axis([0,16,1.e-6,0.1]);
   %loglog(r, abs(u1(:,1)), r, 0.08./r); axis([5,20,1.e-4,0.1]);

%-- 1-D, just in case

   %u0 = 3^(1/4);
   %u  = lsode('udot1d', [u0,0], r);  
   %plot(r, u(:,1), r, u0./sqrt(cosh(2*r)));

%-- number of particles and hamiltonian

   N  = min(find(r>rmax/2)); % truncate data   
   dr = [r(1), r(2:N+1)-r(1:N)];
   r  = [r(1), r(2:N+1)+r(1:N)]/2;
   u  = [u01+u1(1,1),  (u1(2:N+1,1) + u1(1:N,1))']/2;
   v  = [0, (u1(2:N+1,1) + u1(1:N,1))']/2;
   uu = u.^2;
   vv = v.^2;
   
   Ni  = uu .* r;
   KEi = vv.*r;
   PEi = -0.5*uu.^2.*r;
   Pi  = vv .*r.^3;
   Mi  = uu .*r.^3;
   
   NP = 2*pi*sum(Ni.*dr)
   KE = 2*pi*sum(KEi.*dr)
   PE = 2*pi*sum(PEi.*dr)
   P  = 2*pi*sum(Pi.*dr)
   M  = 2*pi*sum(Mi.*dr)/4


%-- print profile to the file

   %drOut=1.e-2;         % density of data in output
   %dN=drOut/dr;

   outfile = fopen(fname, 'wt');
   fprintf(outfile, '#_1.r   2.|psi|   3.d|psi|/dr \n\n');
   fprintf(outfile, '%8.4f  %12.4e  %12.4e\n', 0, u02, 0);

   for i = 2:length(r) % dN *(1:N/dN)
     fprintf(outfile,'%8.4e  %12.4e  %12.4e\n', r(i), u1(i,1), u1(i,2));
   end

   fclose(outfile);

end

%---------------------

function du = udot1d(r,u)    % octave use udot(u,r);  matlab use udot(r,u)

  du = [0 0]';
  
  du(1) =  u(2);
  du(2) =  u(1) - u(1)^5;
 
end

%---------------------

function du = udot2d(u,r)   % octave use udot(u,r);  matlab use udot(r,u)

  b4 = u(3)/4;
  b4=0;
  
  du = [0 0]';
  
  du(1) =  u(2);
  du(2) =  u(1).*(1-b4*r.^2) - u(1)^3 - u(2)/r;
  du(3) =  0;

end

%---------------------

function v = ini(u0,r,beta)

  a  = (1/4)*u0*(1-u0^2);
  b  = (1/16)*a*(1-3*u0^2) - 0.25*beta*u0;

  v1 = u0 + a*r^2 + b*r^4;
  v2 = 2*a*r + 4*b*r^3;

  v  = [v1, v2, beta];
 
end


%---------------------
