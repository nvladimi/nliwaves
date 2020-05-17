function twomode_core(fbase, fnum, Gamma, Rflux, dt, isave, nsave, showplot)
%
% "twomode_core" is a compute core for two-mode evolution, do not modify.
%
% For testing and debugging use "twomode_test".
% For production runs use "twomode_wrap" with "twomode_call" to pass arguments from shell scripts.
%
% Parameters:
%
% fbase        sring base for output files
% seed/fnum    when positive, second argument to "twomode_core" restore IC from (fnum-1) files
%              otherwise use as seed to create (fnum=0) file
% Gamma        strength of decay (<0) or multiplicative forcing (>0), Gamma[re1, im1, re2, im2]
% Rflux        flux of random force,  Rflux[re1, im1, re2, im2]
% dt           timestep
% isave        save data every "isave" timestep
% nsave        generate "nsave" number of saves
% showplot     supress debudding plots if showplot=0
%
  
  
global gamma; global force;

gamma = Gamma;
force = sqrt(Rflux) * sqrt(dt) / sqrt(2);


A = zeros(nsave, 4);
  
%-------------------------------------------------

%-- initial conditions --

   if (fnum > 0) 

       fbase_ic = [fbase, '.',  num2str(fnum-1, '%04d')]

       fbase    = [fbase, '.',  num2str(fnum, '%04d')];

       [f0, t0] = restore_restart(fbase_ic);

       save([fbase, '.param'], 'fnum', 'Gamma', 'Rflux', 'dt', 'isave', 'nsave'); 
  
   else

      save([fbase, '.param'], 'fnum', 'Gamma', 'Rflux', 'dt', 'isave', 'nsave'); 
     
      fbase    = [fbase, '.0000'];
      seed  = -fnum;

      if ( abs(Gamma(1)) ~= 0 )
	   b1 = sqrt( 0.5* Rflux(1)/abs(Gamma(1)) );
      else
	   b1 = 1;
      end

      if ( abs(Gamma(2)) ~= 0 )
	   b2 = -0.5 * sqrt( 0.5* Rflux(2)/abs(Gamma(2)) );
      else
	   b2 = -0.5;
      end


      f0 = [b1, 0, b2, 0];
      t0 = 0;

      randn('twister', seed);

    end
    

	
%-- intergrate ---

  t1 = [0:1]*dt;

   for i2=1:nsave  % cycle over saves

         for i1=1:isave  % cycle over random inputs
   
	    f1 = lsode(['fdotB'; 'fjacB'], f0, t1);
            f2 = f1(end,:);
            f0 = frand(f2, dt);

         end
		 
         A(i2,:) = f0;	 

   end
   
%-- save data --

  t = transpose((1:nsave))*isave*dt + t0;
   
  save_restart(fbase, f0, t(end));

  A = [t,A];
  A = reshape(A, [nsave*5, 1]);

  fid = fopen([fbase, '.a1a2'], 'wb');
  fwrite(fid, A, 'double');
  fclose(fid);


if (showplot == 0) return; end

%-- debugging -- 


   A = reshape(A, [nsave, 5]);
   A = A(:, 2:5);

   [N1, N2, Hnl] = hamiltonian(A);


   figure(1);
   plot(t, A(:,1), '-or', t, A(:,2), '-xr', t, A(:,3), '-ob',  t, A(:,4), '-xb' );
   %axis([t(1), t(end), -2,2]);
   set(gca, "fontsize", 20);

   figure(2);
   N0 = N1(1) + 2*N2(1);
   H0 = Hnl(1);
   plot(t, N1+2*N2 - N0, '-or', t, Hnl - H0, '-ob');
   set(gca, "fontsize", 20)
   %axis([t(1), t(end), -1e-7,1e-7]);
   grid("on");


end



%---------------------

function f = fdotB(x,t)

  global gamma;

  f1 =   2*( x(1)*x(4) - x(2)*x(3) ) + gamma(1)*x(1);
  f2 = - 2*( x(1)*x(3) + x(2)*x(4) ) + gamma(2)*x(2);

  f3 =   2*x(1)*x(2)                 + gamma(3)*x(3);
  f4 = - ( x(1)*x(1) - x(2)*x(2) )   + gamma(4)*x(4);

  f = [f1, f2, f3, f4];

end


%---------------------

function J = fjacB(x,t)

  global gamma;

  J(1,1) =      2 * x(4) + gamma(1);
  J(1,2) =    - 2 * x(3);
  J(1,3) =    - 2 * x(2);
  J(1,4) =      2 * x(1);

  J(2,1) =    - 2 * x(3);
  J(2,2) =    - 2 * x(4) + gamma(2);
  J(2,3) =    - 2 * x(1);
  J(2,4) =    - 2 * x(2);

  J(3,1) =      2 * x(2);
  J(3,2) =      2 * x(1);
  J(3,3) =      gamma(3);
  J(3,4) =  0;

  J(4,1) =    - 2 * x(1);
  J(4,2) =      2 * x(2);
  J(4,3) =  0;
  J(4,4) =      gamma(4);

end


%---------------------

function f = frand(x,dt)

  global force;

  r = randn(1,4) .* force; 

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

   randn('twister');
   randstate = randn ("state");

   fid = fopen([fbase, '.restart'], 'rb');
   randstate = fread(fid, 625, 'uint32');
   f0 = fread(fid, 4, 'double');
   t0 = fread(fid, 1, 'double');
   fclose(fid);

   randn ('state', randstate);

   f0 = reshape(f0,[1,4]);

end

%---------------------

 function [N1, N2, Hnl] = hamiltonian(x);

  N1  = x(:,1).*x(:,1) + x(:,2).*x(:,2);
  N2  = x(:,3).*x(:,3) + x(:,4).*x(:,4) ;
  Hnl  = 2 * ( x(:,1).*x(:,1).*x(:,3) - x(:,2).*x(:,2).*x(:,3) + 2*x(:,1).*x(:,2).*x(:,4) ) ;

 
end


%---------------------
