function fibo_core(fbase, fnum, G, P, m, runtype, dt, isave, nsave, showplot)
%
% "fibo_core" is a compute core for two-mode evolution, do not modify it.
%
% For testing and debugging use "twomode_test".
% For production runs use "twomode_wrap" with "twomode_call" to pass arguments from shell scripts.
%
% Parameters:
%
% fbase        string base for output files
% seed/fnum    when positive, use second argument to "fibo_core" to restore IC from (fnum-1) files,
%              otherwise use it as a seed to create (fnum=0) file
% G            level of damping (application depends on runtype)
% P            level of forcing (application depends on runtype)
% m            number of modes
% runtype      integer to distinguish different forcing and damping and IC
% dt           timestep
% isave        save data every "isave" timestep
% nsave        generate "nsave" number of saves
% showplot     supress debudding plots if showplot=0
%
  
  
global gamma;      % gamma(1,M)    complex array
global force;      % force(1,2*M)  real array 
global M;          % number of modes, M=m

M = m;

%-------------------------------------------------


%-- initial conditions --

   if (fnum > 0) 

      fbase_ic = [fbase, '.',  num2str(fnum-1, '%04d')];

      fbase    = [fbase, '.',  num2str(fnum, '%04d')];

      [f0, t0] = restore_restart(fbase_ic);

      save([fbase, '.param'], 'fnum', 'G', 'P', 'm', 'runtype', 'dt', 'isave', 'nsave'); 
  
   else

      save([fbase, '.param'], 'fnum', 'G', 'P', 'm', 'runtype', 'dt', 'isave', 'nsave'); 
     
      seed  = -fnum;
      fnum  = 0;
      fbase    = [fbase, '.0000'];
      save([fbase, '.param'], 'fnum', 'G', 'P', 'm', 'runtype', 'dt', 'isave', 'nsave'); 

      randn('twister', seed);
      rand('twister', seed);
 
      f0 = set_IC(G, P, runtype);

      t0 = 0;

    end

    set_forcing_damping(G, P, runtype);

    force = force * sqrt(dt);

    A = zeros(nsave, 2*M);

	
%-- intergrate ---

  t1 = [0:1]*dt;

   for i2=1:nsave  % cycle over saves

         for i1=1:isave  % cycle over random inputs
   
  	    f1 = lsode(['fdotA'; 'fjacA'], f0, t1);
            f2 = f1(end,:);
            f0 = frand(f2, dt);

         end
		 
         A(i2,:) = f0;	 

   end
   
%-- save data --

  t = transpose((1:nsave))*isave*dt + t0;
   
  save_restart(fbase, f0, t(end));

  A = [t,A];
  A = reshape(A, [nsave*(2*M+1), 1]);

  fid = fopen([fbase, '.ak'], 'wb');
  fwrite(fid, A, 'double');
  fclose(fid);


if (showplot == 0) return; end

%-- debugging -- 


   A = reshape(A, [nsave, 2*M+1]);
   ind = (1:M)*2;
   reA = squeeze(A(:,ind));
   imA = squeeze(A(:,ind+1));
   absA = sqrt(reA.*reA + imA.*imA);

   t=t/pi;

   figure(1);
   plot(t, reA(:,1), '-+r', t, imA(:,2), '-r', ...
        t, reA(:,2), '-xb', t, imA(:,2), '-b', ...
        t, reA(:,3), '-ok', t, imA(:,3), '-k');
   set(gca, "fontsize", 20); grid("on");

   figure(2);
   plot(t, absA, '-k', t, absA(:,1), '+r',  t, absA(:,2), 'xb');
   set(gca, "fontsize", 20); grid("on");

end

%---------------------

function f = fdotA(x,t) 
  
  global gamma;
  global M;

  q=[[0;0;0;0];x;[0;0;0;0]];
  q = reshape(q,2,M+4);

  a = q(1,:) + 1i*q(2,:);

  am2 = a(1:M);
  am1 = a(2:M+1);
  a0  = a(3:M+2);
  ap1 = a(4:M+3);
  ap2 = a(5:M+4);


  q = -1i * ( am2 .* am1  + conj(am1) .* ap1 + conj(ap1) .* ap2) + gamma.*a0; 

  f = [real(q); imag(q)];
  f = reshape(f, 1, 2*M);

end


%---------------------

function J = fjacA(x,t)

  global gamma;
  global M;
 
  q =[[0;0;0;0];x;[0;0;0;0]];
  q = reshape(q,2,M+4);
  x = q(1,:);
  y = q(2,:);

  J = zeros(2*M, 2*M+8);

  for m=1:M
     k=m+2;
     ind = 2*m-1:2*m+8;
     u = [ y(k-1), x(k-1),  y(k-2)+y(k+1), x(k-2)-x(k+1), 0, 0, -y(k-1)+y(k+2),  x(k-1)-x(k+2), -y(k+1),  x(k+1)]; 
     v = [-x(k-1), y(k-1), -x(k-2)-x(k+1), y(k-2)-y(k+1), 0, 0, -x(k-1)-x(k+2), -y(k-1)-y(k+2), -x(k+1), -y(k+1)]; 
     J(2*m-1,ind) = u;
     J(2*m,  ind) = v;
  end
  
  J = J(:, 5:2*M+4);

end

%---------------------
%---------------------

function f = frand(x,dt)

  global force;
  global M;

  r = randn(1, 2*M) .* force; 

  f = x + r;

end


%---------------------

function save_restart(fbase, f0, t0)

   randstate = randn ("state");

   fid = fopen([fbase, '.restart'], 'wb');
   fwrite(fid, randstate, 'uint32');
   fwrite(fid, f0, 'double');
   fwrite(fid, t0, 'double');
   fclose(fid);

end

%---------------------


function [f0, t0] = restore_restart(fbase)

   global M;
  
   randn('twister');
   randstate = randn ("state");

   load([fbase, '.param']);
  

   fid = fopen([fbase, '.restart'], 'rb');
   randstate = fread(fid, 625, 'uint32');
   f0 = fread(fid, 2*m, 'double');
   t0 = fread(fid,   1, 'double');
   fclose(fid);

   randn ('state', randstate);

   f0 = reshape(f0,[1, 2*m]);

   if (m < M)
     disp('extending number of modes');
     f0 = [f0, zeros(1, 2*(M-m))];
   end

   if (m > M)
     disp('reducing number of modes');
     f0 = f0(1:2*M);
   end
     
end

%---------------------

function  f0 = set_IC(G, P, runtype)
%
%  first digit in runtype determines initial condition
%  second digit determines forcing
%
%  ic_type 1:     a0 =  [1, 1, 0, ..., 0]
%  ic_type 2:     a0 =  [1, 1, -1/sqrt(2), 0, ..., 0]
%  ic_type 3:     |a0|^2 = 1/Fi,   random phases
%  ic_type 4:     |a0|^2 = random, random phases 
%
     
   global M;

   Fi = fibonacci(M);

   ic_type = fix(runtype/10);

   if (ic_type == 1)
     
     f0 = zeros(2,M);
     f0(1,1) = 1;
     f0(1,2) = 1;
 
   elseif (ic_type == 2)

     f0 = zeros(2,M);
     f0(1,1) = 1;
     f0(1,2) = 1;
     f0(1,3) = -1/sqrt(2);

   elseif (ic_type == 3)

     amp = 1./sqrt(Fi);
     phi = rand(1,M)*2*pi;

     f0(1,:) = amp.*cos(phi);
     f0(2,:) = amp.*sin(phi);

   elseif (ic_type == 4)

     amp = sqrt(rand(1,M));
     phi = rand(1,M)*2*pi;
     f0(1,:) = amp.*cos(phi);
     f0(2,:) = amp.*sin(phi);

   else 
     disp('Unknown type of initial conditions');
     f0 = zeros(2,M);
   end

   f0 = reshape(f0, [1,2*M]);

   % f0 = [sqrt(2), 0, sqrt(2), 0, -1, 0];
   % f0 = 2*rand(1,2*M) - 1;

end
%---------------------


function   set_forcing_damping(G, P, runtype)
%
%  first digit in runtype determines initial condition
%  second digit determines forcing
%
%  force_type 0:     force = 0, gamma = 0;
%  force_type 1:     force = [P, P, 0, ..., 0],         gamma = [0, ..., 0, G, G];        direct cascade
%  force_type 2:     force = [P, P, 0, ..., 0],         gamma = [0, ..., 0, G*Fi, G*Fi];  direct cascade balanced  
%  force_type 3:     force = [0, ..., 0, P, P],         gamma = [G, G, 0, ..., 0];        inverse cascade
%  force_type 4:     force = [0, ..., 0, P/Fi, P/Fi],   gamma = [G, G, 0, ..., 0];        inverse cascage balanced
%  force_type 5:     force = [P, P, P, ..., P],         gamma = P*Fi/2;                   equilibrium by force
%  force_type 6:     force = 2*G/Fi,                    gamma = [G, G, G, ..., G];        equilibrium by decay
%

   global M;
   global gamma;      % gamma(1,M)    complex array
   global force;      % force(1,2*M)  real array 

   Fi = fibonacci(M);

   gamma = zeros(1,M); 
   force = zeros(1,M);

   force_type = rem(runtype,10);

   if (force_type == 0)
     
     % do nothing;

   elseif (force_type == 1)  % direct cascade

     force(1) = P;
     force(2) = P;
     gamma(end) = G;
     gamma(end-1) = G;

   elseif (force_type == 2)  % balanced direct cascade
				  
     force(1) = P;
     force(2) = P;
     gamma(end) = G*Fi(end);
     gamma(end-1) = G*Fi(end-1);

   elseif (force_type == 3)  % inverse cascade

     force(end) = P;
     force(end-1) = P;
     gamma(1) = G;
     gamma(2) = G;

   elseif (force_type == 4)  % balanced inverse cascade

     force(end) = P/Fi(end);
     force(end-1) = P/Fi(end-1);
     gamma(1) = G;
     gamma(2) = G;

   elseif (force_type == 5)  % constant force

     force = ones(1,M)*P;
     gamma = 0.5 * force .* Fi;

   elseif (force_type == 6)  % constant damping

     gamma = ones(1,M)*G;
     force = 2 * gamma ./ Fi;

   else 
     disp('Unknown type of forcing');
   end
  
   gamma = -gamma;
   force = sqrt([force; force]) / sqrt(2);
   force = reshape(force, [1,2*M]);

end

%---------------------
