%function [x,u_out] = yang_ncg  % manual selection of cutoff
function yang_ncg
%
% modified from
% Jianke Yang, p8.m: the Newton-CG method for computing gap solitons 
% in the 2D NLS equation with a potential and self-focusing nonlinearity:
%
% u_{xx}+u_{yy} - V(x,y) u + u^3 = lambda*u  
%
% V(x,y) = -beta/4 (x^2 + y^2)

global rr
global kk
global Mv
global u
global lambda
global mu
global nu


fn0= -5;    r0=36.0;  u0=2.206;                  fbase='mu0_a1';   follow_tail = false;


%follow_tail = true;
%
%fn0= 295;  r0=36.7;  u0=2.190;  drdbeta= -350;  fbase='mu0_a2';
%fn0= 495;  r0=35.8;  u0=2.179;  drdbeta= -180;  fbase='mu0_a3';
%fn0= 795;  r0=36.1;  u0=2.160;  drdbeta= -110;  fbase='mu0_a4';
%fn0=1295;  r0=36.2;  u0=2.124;  drdbeta=  -65;  fbase='mu0_a5';   % also "b5", at N=640 

 fnstep =  5;
 fnmax  =  2000;
 dbeta  =  0.0001 * fnstep;
 beta   =  0.0001 * fn0;
 r      =  r0;

 Lx =  80; 
 Ly =  80; 
 N  = 320; 

%----------------------------------------------
 
 lambda = 1;
 mu =     0;
 nu =     0.5;


%----------------------------------------------


 x=-Lx/2:Lx/N:Lx/2-Lx/N; kx=[0:N/2-1  -N/2:-1]*2*pi/Lx;
 y=-Ly/2:Ly/N:Ly/2-Ly/N; ky=[0:N/2-1  -N/2:-1]*2*pi/Ly;

 [X,Y]=meshgrid(x,y); [KX,KY]=meshgrid(kx,ky);

 kk = KX.^2+KY.^2;
 rr = X.^2 + Y.^2;
 Mv = KX.*KX ./ (KX.*KX + nu*KY.*KY);  Mv(1,1) = 0;

 clear KX; clear KY;

 fid = fopen([fbase, '.dat'], 'at');
 fprintf(fid, '#1.beta  2.numwaves  3.h  4.phi  5.tail_r  6.tail_u  7.err  8.newton_steps  9.CG_iters\n\n');
 fclose(fid);


 %-- manual selection of cutoff --
 %
 %  u = u0*exp(-rr/2);  for i=1:4  allout=Newton_iter(r0, beta); end;  u_out= u(end/2+1,:);  return


%--------------------

 for fn=fn0+fnstep:fnstep:fnmax

    beta  =  beta + dbeta;

    allall=[];


if follow_tail 

    % Dr = pi / (sqrt(beta)*r0/4 - 1/sqrt(beta)*log(r0)/r0);

    %-- find minimum --


    dr =  drdbeta * dbeta;
    rrange = r0 +  (0:2) * dr;

    for r = rrange

       u = u0*exp(-rr/2);

       allout = Newton_iter(r, beta);   iters     = allout(1);
       allout = Newton_iter(r, beta);   iters     = iters + allout(1);
       allout = Newton_iter(r, beta);   iters     = iters + allout(1);
       allout = Newton_iter(r, beta);   allout(1) = allout(1) + iters;
       allall = [allall; allout];
       
    end


   allu = allall(:,3);
   allr = allall(:,2);

   p = polyfit (allr, allu, 2);

   r = - 0.5*p(2)/p(1);

end 

   %-- final iteratons --
 
   u = u0*exp(-rr/2);

   allout = Newton_iter(r, beta);   allall = [allall; allout];
   allout = Newton_iter(r, beta);   allall = [allall; allout];
   allout = Newton_iter(r, beta);   allall = [allall; allout];
   allout = Newton_iter(r, beta);   allall = [allall; allout];
 
   utail = allout(3);
   u0    = allout(5);

   num_waves = sum(sum(u.*u)) *Lx*Ly/N/N;

   %figure(1); plot(allr, allu, r, utail0, '*', r, utail, 'o') ; grid("on"); 


   %-- record data --

    fid = fopen([fbase, num2str(fn,'.%04d'), '.dat'], 'wt');
       fprintf(fid, '#1.step  2.cg_iters   3.tail_r  4.tail_u  5.err  6.u_center  7.phi_center\n\n');
       for i=1:length(allall(:,1))
         fprintf( fid, ' %3d   %6d     %10.6f   %10.6f   %10.4e  %10.6f  %10.6f\n', ...
		  i,   allall(i,:) );
       end
    fclose(fid);


    fid = fopen([fbase, '.dat'], 'at');
       fprintf(fid, ' %6.4f   %10.6f  %14.10f  %10.6f  %10.6f  %10.6f  %10.4e  %4d  %6d\n', ...
	       beta,  num_waves,  ...
               allall(end,5),  allall(end,6),  allall(end,2),  allall(end,3),  allall(end,4), ...
               length(allall(:,1)),  sum(allall(:,1)) );
    fclose(fid);


  %  fid = fopen([fbase, num2str(fn,'.%04d')] , 'wb');
  %      fwrite(fid, u, 'double');
  %  fclose(fid);


  %--  new guess --

  drdbeta = (r - r0) / dbeta;
  r0 = r;

  end

  %-- loop over beta ends 

end



%-----------------------------------

function allout = Newton_iter(r0, beta)

global rr
global kk
global Mv
global u
global lambda
global mu
global nu

       errorCG  = 1e-2;
       c = 3;

       %---------------------

       pow = 60;
       V = - 0.25 * beta * rr .* exp (-(sqrt(rr)/r0).^pow);

       phi = ifft2(  Mv.*fft2(u.*u) );
       L0u = ifft2( -kk.*fft2(u)) - (V - u.*u + lambda).*u - mu*u.*phi;
        
       %-- CG iterations begin

       ncg   = 0;
       du    = zeros(size(kk));
       R     = -L0u; 
       MinvR = ifft2( fft2(R)./(c+kk) );
       R2    = sum(sum(R.*MinvR)); 
       R20   = R2;
       Du    = MinvR;
 
       while ( R2 > R20*errorCG^2 )

          dphi  = ifft2( Mv.*fft2(2.*u.*Du) );
      
 	  L1u   = ifft2(-kk.*fft2(Du)) - (V - 3*u.*u + lambda).*Du ...
	             - mu*phi.*Du - mu*u.*dphi;

          a     = R2/sum(sum(Du.*L1u));
          du    = du + a*Du;
          R     = R  - a*L1u; 
          MinvR = ifft2( fft2(R)./(c+kk) );
          R2old = R2;
          R2    = sum(sum(R.*MinvR));
          b     = R2/R2old;
          Du    = MinvR + b*Du;
          ncg   = ncg + 1;

       end

       %-- CG iterations end --

       u = u + real(du);

       u0    = u(end/2+1,:);
       x0    = sqrt(rr(end/2+1,:));
       du0   = abs(u0)-circshift(abs(u0),[0,1]);
       ind   = max(find(du0>1.e-6));
       utail = abs(u0(ind));

       % plot(x0,u0, '.-', x0(ind), u0(ind), 'or'); axis([20,40,-0.003, 0.003]); grid('on')


       phi = ifft2(  Mv.*fft2(u.*u) );
       L0u = ifft2( -kk.*fft2(u)) - (V - u.*u + lambda).*u - mu*u.*phi;

       err = max(max(abs(L0u) ./ (1+rr) )); 

       ucenter =  u(end/2+1, end/2+1);
       phicenter = real(phi(end/2+1, end/2+1));

       allout = [ncg, r0, utail, err, ucenter, phicenter];
 
end

%----------------------------------------------------------------
