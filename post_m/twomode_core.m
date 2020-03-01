function twomode_core(fbase, fnum, Epsilon, Gamma, Theta, n)

global epsilon; global gamma; global theta;

epsilon = Epsilon; gamma = Gamma;  theta = Theta;

n1 = n(1);         % number of timesteps per t=pi, must be multiple of n2
n2 = n(2);         % number of random inputs per t=pi, must be multiple of n3
n3 = n(3);         % number of data saves per t=pi
n4 = n(4);         % number of pi-periods to run
N  = n3*n4;        % total number of snapshots to save

A = zeros(N, 4);

  
%-------------------------------------------------

%-- initial conditions --

   if (fnum > 0) 

       fbase_ic = [fbase, '.',  num2str(fnum-1, '%04d')]

       fbase    = [fbase, '.',  num2str(fnum, '%04d')];

       [f0, t0] = restore_restart(fbase_ic);

       save([fbase, '.param'], 'fbase', 'fnum', 'epsilon', 'gamma', 'theta', 'n'); 
  
   else

      save([fbase, '.param'], 'fbase', 'fnum', 'epsilon', 'gamma', 'theta', 'n'); 
     
      fbase    = [fbase, '.0000'];
      seed  = -fnum;
     
%%%   b1 = 2/ epsilon;   b2 = - 1/epsilon;        %-- including linear terms --
      b1 = 1/ epsilon;   b2 = - 0.5/epsilon;      %-- excluding linear terms --

      f0 = [b1, 0, b2, 0];
      t0 = 0;

      rand('twister', seed);

    end
    

	
%-- intergrate ---

   dt1 = pi/n1;
   dt2 = pi/n2;


   t1 = (0:n1/n2)*dt1;

   for i3=1:N  % cycle over saves

         for i1=1:n2/n3  % cycle over random inputs
   
%	    f1 = lsode('fdotB', f0, t1);
	    f1 = lsode(['fdotB'; 'fjacB'], f0, t1);
            f2 = f1(end,:);
            f0 = frand(f2, dt2);

         end
		 
         A(i3,:) = f0;	 

   end
   
%-- save data --

  dt3 = pi/n3;
  t = transpose((1:N))*dt3 + t0;
   
  save_restart(fbase, f0, t(end));

  A = [t,A];
  A = reshape(A, [N*5, 1]);

  fid = fopen([fbase, '.a1a2'], 'wb');
  fwrite(fid, A, 'double');
  fclose(fid);


%return

%-- debugging -- 


   A = reshape(A, [N, 5]);
   A = A(:, 2:5);

   [N1, N2, Hnl] = hamiltonian(A);


   figure(1);
   plot(t, A(:,1), '-or', t, A(:,2), '-xr', t, A(:,3), '-ob',  t, A(:,4), '-xb' );
   axis([t(1), t(end), -20,20]);
   set(gca, "fontsize", 20);

   figure(2);
   plot(t, N1+2*N2 - 1.5*b1*b1, '-or', t, Hnl + b1*b1, '-ob');
   set(gca, "fontsize", 20)
   axis([t(1), t(end), -1e-5,1e-5]);
   grid("on");


end

%---------------------

function f = fdotA(x,t)

  global epsilon;
  global gamma;

  f1 =     x(2) + 2*epsilon*( x(1)*x(4) - x(2)*x(3) ) + gamma(1)*x(1);
  f2 =   - x(1) - 2*epsilon*( x(1)*x(3) + x(2)*x(4) ) + gamma(2)*x(2);

  f3 =   2*x(4) + 2*epsilon* x(1)*x(2)                + gamma(3)*x(3);
  f4 = - 2*x(3) -   epsilon*( x(1)*x(1) - x(2)*x(2) ) + gamma(4)*x(4);

  f = [f1, f2, f3, f4];

end


%---------------------

function f = fdotB(x,t)

  global epsilon;
  global gamma;

  f1 =   2*epsilon*( x(1)*x(4) - x(2)*x(3) ) + gamma(1)*x(1);
  f2 = - 2*epsilon*( x(1)*x(3) + x(2)*x(4) ) + gamma(2)*x(2);

  f3 =   2*epsilon* x(1)*x(2)                + gamma(3)*x(3);
  f4 =  - epsilon*( x(1)*x(1) - x(2)*x(2) )  + gamma(4)*x(4);

  f = [f1, f2, f3, f4];

end


%---------------------

function J = fjacA(x,t)

  global epsilon;
  global gamma;

  eps2 = 2*epsilon;

  J(1,1) =      eps2 * x(4) + gamma(1);
  J(1,2) =  1 - eps2 * x(3);
  J(1,3) =    - eps2 * x(2);
  J(1,4) =      eps2 * x(1);

  J(2,1) = -1 - eps2 * x(3);
  J(2,2) =    - eps2 * x(4) + gamma(2);
  J(2,3) =    - eps2 * x(1);
  J(2,4) =    - eps2 * x(2);

  J(3,1) =      eps2 * x(2);
  J(3,2) =      eps2 * x(1);
  J(3,3) =      gamma(3);
  J(3,4) =  2;

  J(4,1) =    - eps2 * x(1);
  J(4,2) =      eps2 * x(2);
  J(4,3) =  -2;
  J(4,4) =      gamma(4);

end



%---------------------

function J = fjacB(x,t)

  global epsilon;
  global gamma;

  eps2 = 2*epsilon;

  J(1,1) =      eps2 * x(4) + gamma(1);
  J(1,2) =    - eps2 * x(3);
  J(1,3) =    - eps2 * x(2);
  J(1,4) =      eps2 * x(1);

  J(2,1) =    - eps2 * x(3);
  J(2,2) =    - eps2 * x(4) + gamma(2);
  J(2,3) =    - eps2 * x(1);
  J(2,4) =    - eps2 * x(2);

  J(3,1) =      eps2 * x(2);
  J(3,2) =      eps2 * x(1);
  J(3,3) =      gamma(3);
  J(3,4) =  0;

  J(4,1) =    - eps2 * x(1);
  J(4,2) =      eps2 * x(2);
  J(4,3) =  0;
  J(4,4) =      gamma(4);

end


%---------------------

function f = frand(x,dt)

  global theta;

  r = randn(1,4) .* theta; 

  f = x + r * sqrt(dt);

end


%---------------------

function save_restart(fbase, f0, t0)

   randstate = rand ("state");

   fid = fopen([fbase, '.restart'], 'wb');
   fwrite(fid, randstate, 'uint32');
   fwrite(fid, f0, 'double');
   fwrite(fid, t0,  'double');
   fclose(fid);

end

%---------------------


function [f0, t0] = restore_restart(fbase)

   rand('twister');
   randstate = rand ("state");

   fid = fopen([fbase, '.restart'], 'rb');
   randstate = fread(fid, 625, 'uint32');
   f0 = fread(fid, 4, 'double');
   t0 = fread(fid, 1, 'double');
   fclose(fid);

   rand ('state', randstate);

   f0 = reshape(f0,[1,4]);

end

%---------------------

 function [N1, N2, Hnl] = hamiltonian(x);

  global epsilon;

  N1  = x(:,1).*x(:,1) + x(:,2).*x(:,2);
  N2  = x(:,3).*x(:,3) + x(:,4).*x(:,4) ;
  Hnl  = epsilon * 2 * ( x(:,1).*x(:,1).*x(:,3) - x(:,2).*x(:,2).*x(:,3) + 2*x(:,1).*x(:,2).*x(:,4) ) ;

 
end


%---------------------
