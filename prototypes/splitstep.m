function splitstep
% solves NLS equation 
%
%   i dPsi/dt + (1-i*a*eps)*del^2(Psi) + (1+i*eps)*|Psi|^2*Psi = i*b*eps*Psi 
%
% using Dyachenko-Newell-Pushkarev-Zakharov-PhysicaD-1992 algorithm

%-- Input parameters --  

  a       =  1;             % \ 
  b       =  20;             % -  coefficients in the equation 
  epsilon =  0.1;           % /

  L     =   12.8;           % length of the box
  r0    =    1.0;           % initial width of gaussian psi
  h0    =    4.0;           % initial height of gaussian psi
  x0    =    0.0;           % initial x-location of gaussian psi
  y0    =    0.0;           % initial y-location of gaussian psi      

  tmax  =    0.1;           % maximum time
  tinc  =    0.01;          % how often to post-process

  N     =    256;           % discretization points 

  dtmax =    0.001;          % maximum timestep

  fbase = 'OUT/dnpz1024_dt1e-3_da';       % name base for output files
  fbase = 'OUT/test';

%-- Initializaion --
 
  %[U, x] = initPsi_1D(L, N, h0, r0);
  [U, x]  = initPsi_2D(L, N, h0, r0, x0, y0);

  fU = fftn(U);
  fU = dealias(fU);
 
  DD = laplacian(N,L,dims(U));  % matrix of wavenumbers squared
  Q  = (a*epsilon + i)*DD  + b*epsilon;

  t     = 0;
  iter  = 0;
  tplot = tinc;
  dx    = L/N;

  if dims(U) == 1
    postproc_1d(U,x,t,fbase);
  else
    postproc_2d(U,x,t,fbase);
  end

%-- Main loop - time advancement --

  while t<tmax

    dt = min(1/max(max(abs(U).^2)), dtmax);
    if (t+dt > tplot) 
        dt = tplot-t;
    end 
    dt = dtmax;
  
    %------------------------

    % advance linear part half-step
    fU = (1 + 0.25*dt*Q)./(1 - 0.25*dt*Q) .* fU; 
    U  =  ifftn(fU);
 
    % advance non-linear part whole step

    UU = (abs(U)).^2;
    aa =  (i-epsilon)/(2*epsilon);
    U  = U .* exp( aa * log(1 + 2*epsilon*UU*dt) );

    % advance linear part half-step
    fU = fftn(U);
    %fU = dealias(fU);
    fU = (1 + 0.25*dt*Q)./(1 - 0.25*dt*Q) .* fU;
 
    t = t + dt;

    %------------------------


    % plot and some diagnostics

    if abs(t-tplot)<eps

      U  =  ifftn(fU);

      %fprintf('t= %f,  dt= %f,  s= %d,  iter=%d\n', t, dt, s, iter);

      if dims(U)==1
        postproc_1d(U,x,t,fbase);
        plot(x, real(U), x, imag(U), x, abs(U));  axis([-L/2,L/2,-4,4]);
        pause(0.1);
      else
        postproc_2d(U,x,t,fbase);
        %mesh(x, x, abs(U));       
        %pause(0.1);
      end
    
      tplot = tplot + tinc;

      U  =  fftn(fU);

    end
    
  end

  % show final

  if dims(U)==2
    %mesh(x, x, abs(U));
    %postproc_2d(U,x,t,fbase);
  end

end

%-------------------------------------------------
function [U,x] = initPsi_1D(L, N, h0, r0)

  dx = L/N;

  x  = (-N/2:1:N/2-1)*dx + dx/2;
  U  = h0*exp( -(x/r0).^2 ) + 0*i;
  %plot(x, real(U), x, imag(U), x, abs(U));

end

%-------------------------------------------------
function [U,x] = initPsi_1G(L, N, h0, r0)

  dx = L/N;

  x  = (-N/2:1:N/2-1)*dx + dx/2;
  U  = h0 ./ cosh(x*h0/sqrt(2)) + 0*i;
  %plot(x, real(U), x, imag(U), x, abs(U));

end

%-------------------------------------------------
function [U,x,y] = initPsi_2D(L, N, h0, r0, x0, y0)

  dx = L/N;

  x  = (-N/2:1:N/2-1)*dx + dx/2;
  y  = (-N/2:1:N/2-1)*dx + dx/2;

  [xx, yy] = meshgrid (x, y);

  rr = (xx-x0).^2 + (yy-y0).^2;

  U  = h0*exp( -rr/r0^2 ) + 0*i;
  %mesh(x, y, real(U));

end

%-------------------------------------------------

function DD = laplacian(N,L,d)
% create a matrix of wavenumbers squared 
% normalized to the domain size

% create the vector of k 
% k = [0,1,2...,N/2, -N/2, ..., -2, -1]
   
  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);

  if d==1
    DD = -(k*2*pi/L).^2;
  else
    [k1, k2] = meshgrid (k, k);
    DD = -(k1.^2 + k2.^2)*(2*pi/L)^2;
  end

end

%-------------------------------------------------

function fU=dealias(fUin)
% zero out high harmonics

  fU = fUin;

  N = length(fU);
  d  = dims(fU);

  if d==1
    fU(N/4:3*N/4+1) = 0;
  else
    fU(N/4:3*N/4+1, :) = 0;
    fU(:, N/4:3*N/4+1) = 0;
  end

end

%-------------------------------------------------

function d = dims(U)
% dimensionality of the matrix

  [r, c] = size(U);

  if r==1 
      d=1;
  else
      d=2;
  end

end

%-------------------------------------------------
