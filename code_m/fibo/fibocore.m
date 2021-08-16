function fibocore(runtype, fbase, fnum, alpha, m, mp, PI, G1, G2, dt, isave, nsave, showplot)
%
% Parameters:
%
% runtype      string: "fibA", "fibB", "dblA", or "dblB" 
% fbase        string base for output files
% seed/fnum    when positive, use third argument to "fibo_core" to restore IC from 
%              (fnum-1) file, otherwise use it as a seed to create (fnum=0) file
% alpha        parameter controlling interactions, as in papers
% m            number of modes
% mp           pumped mode
% PI           influx at the pumped mode
% G1           rate of damping at low modes
% G2           rate of damping at high modes
% dt           timestep
% isave        save data every "isave" timestep
% nsave        generate "nsave" number of saves
% showplot     supress debudding plots if showplot=0
%

global RUNTYPE     % string: "fibA", "fibB", "dblA", or "dblB" 
global gamma;      % gamma(1,M)    complex array
global force;      % force(1,2*M)  real array 
global M;          % number of modes, M=m
global V0;         % weights for mode interactions at i
global V1;         % weights for mode interactions at i-1
global V2;         % weights for mode interactions at i-2
  
M = m;
RUNTYPE = runtype;

%-------------------------------------------------


%-- initial conditions --

   if (fnum > 0)

      fbase_ic = [fbase, '.',  num2str(fnum-1, '%04d')];

      fbase    = [fbase, '.',  num2str(fnum, '%04d')];

      [f0, t0] = restore_restart(fbase_ic);

      save([fbase, '.param'], 'fnum', "PI", 'G1', 'G2', 'm', 'alpha', 'runtype', 'dt', 'isave', 'nsave'); 
  
   else 

      save([fbase, '.param'], 'fnum', "PI", 'G1', 'G2', 'm', 'alpha', 'runtype', 'dt', 'isave', 'nsave'); 
     
      seed  = -fnum;
      fnum  = 0;
      fbase    = [fbase, '.0000'];
      save([fbase, '.param'], 'fnum', "PI", 'G1', 'G2', 'm', 'alpha', 'runtype', 'dt', 'isave', 'nsave'); 

      randn('twister', seed);
      rand('twister', seed);
 
      f0 = set_IC(runtype);

      t0 = 0;

    end

%-- interactions --

    if (RUNTYPE == "fibA")

       Fi = fibonacci(M);
       V0  = Fi.^alpha;
       V1 =    [0, V0(1:M-1)];
       V2 = [0, 0, V0(1:M-2)];

       Pp  = PI/Fi(mp);
  
    elseif (RUNTYPE == "fibB")

       i = 1:M;
       phi = (1 + sqrt(5))/2;
       q = i*(2*alpha - 1) + 1 - alpha;
       V1  = (PI * phi.^q).^(1/3)  ;
       V0 =  V1/phi;
       V2 =  V1*phi;

       Pp  = PI;

    elseif (RUNTYPE == "dblA")

       i = 1:M;
       V1  = 2.^(alpha*i - alpha) ;
       V0 =  2.^(alpha*i + 1);

       Pp  = PI/2^mp;

    elseif (RUNTYPE == "dblB")

       i = 1:M;
       q = i*(2*alpha - 1) - alpha;
       V1 = (PI * 2.^q).^(1/3);
       V0 =  V1;

       Pp  = PI;
    
   else
     disp('Unknown runtype');
     return;
   end
    
 
    
%-- forcing and damping --

    set_forcing_damping(mp, Pp*dt, G1, G2);


%-- empty array to store evolution --

    A = zeros(nsave, 2*M);


%-- intergrate ---

  t1 = [0:1]*dt;

   for i2=1:nsave  % cycle over saves

         for i1=1:isave  % cycle over random inputs
   
  	    f1 = lsode('fdot', f0, t1);
            f2 = f1(end,:);
            f0 = f2 + random_forcing();

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

function f = fdot(x,t) 

   global RUNTYPE;
   global M;
   global gamma;
   global V0;
   global V1;
   global V2;

   %-- fold real array into complex --

   q = [[0;0;0;0];x;[0;0;0;0]];
   q = reshape(q,2,M+4);

   a = q(1,:) + 1i*q(2,:);

   %-- arrays a_{i-2},  a_{i-1},  a_{i},  a_{i+1},  a_{i+2} --

   am2 = a(1:M);
   am1 = a(2:M+1);
   a0  = a(3:M+2);
   ap1 = a(4:M+3);
   ap2 = a(5:M+4);


   %-- compute RHS for model-specific ODEs without forcing --

   if ((RUNTYPE == "fibA") | (RUNTYPE == "fibB"))

        q = -1i * ( V2.*am2 .* am1  + V1.*conj(am1) .* ap1 + V0.*conj(ap1) .* ap2) + gamma.*a0; 

   elseif ((RUNTYPE == "dblA") | (RUNTYPE == "dblB"))

        q = -1i * ( V1.* am1 .* am1 + V0.*conj(a0) .* ap1) + gamma.*a0; 
 
   else
     disp('Unknown runtype');
     return;
   end


   %-- unfold comlex array into real --

   f = [real(q); imag(q)];
   f = reshape(f, 1, 2*M);

end


%---------------------

function J = fjac(x,t)
%
% unused jacobian for original "fibo" ODE;
% supposed to make integrations faster and more accurate;
% verify for "fibo" and adjust for "gold", "dlba", and "dblb";
% use as:  lsode(['fdot'; 'fjac'], f0, t1);

  
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

function r = random_forcing()

  global force;
  global M;

  r = randn(1, 2*M) .* force; 

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

function  f0 = set_IC(runtype)
%
%  random phases; amplidudes according to
%  "fibA":     |a0|^2 = 1/Fi
%  "fibB":     |a0|^2 = random [0:1]
%  "dblA":     |a0|^2 = 2^{-i} 
%  "dblB":     |a0|^2 = random [0:1]

   global M;

   if (runtype == "fibA")

     Fi = fibonacci(M);
     amp = 1./sqrt(Fi);

   elseif (runtype == "fibB")

     amp = sqrt(rand(1,M));

   elseif (runtype == "dblA")

     i = 0:M-1; 
     amp = 1./sqrt(2.^i);

   elseif (runtype == "dblB")

     amp = sqrt(rand(1,M));
 
   else
     disp('Unknown type of initial conditions');
     amp = zeros(2,M);
   end

   phi = rand(1,M)*2*pi;
   f0(1,:) = amp.*cos(phi);
   f0(2,:) = amp.*sin(phi);

   f0 = reshape(f0, [1,2*M]);

end
%---------------------


function   set_forcing_damping(mp, P, G1, G2)
%
%  first digit in runtype determines initial condition
%  second two digits determine forcing mode
%

   global M;
   global gamma;      % gamma(1,M)    complex array
   global force;      % force(1,2*M)  real array 


   gamma = zeros(1,M); 
   force = zeros(1,M);

   force(mp)        =  P;
   gamma(1:2)       = -G1;
   gamma(end-1:end) = -G2;
  
   force = sqrt([force; force]) / sqrt(2);
   force = reshape(force, [1,2*M]);

end

%---------------------
function fn = fibonacci(m)

    fn = zeros(1, m);
    fn(1) = 1;
    fn(2) = 1;
    for k=3:m
        fn(k) = fn(k-1)+fn(k-2);
    end
    
end

%---------------------
