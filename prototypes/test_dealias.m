function test_dealias
%
%  Illustration of the effect of dealiasing.
%
%  Solves a simple diffusion equation with U^2 nonlinearity
%
%    dU/dT = Del(U) + R(U)*sign(dU/dx),    R(U) = U*(1-U)/4;
%
%  To see the effect of dealiasing, at the end of this file 
%  uncomment zeroing out of high frequency modes.
%
%  Screen output: plots of U, abs(fU), dU/dx, and d2U/dx2.


%-- Input parameters --  

  N     =   200;            % discretization points 
  L     =    200;           % length of the box
  r0    =   10.0;           % initial width of gaussian psi
  h0    =    1.0;           % initial height of gaussian psi

  tmax  =  100.1;           % maximum time
  dt    =    0.1;           % timestep
  dn    =    100;           % how often to post-process

%-- Initializaion --
 
  %[U, x] = init_gauss(L, N, h0, r0);
  [U, x] = init_tanh(L, N, h0, r0);
  %plot(x, U, '-*');  axis([-L/2,L/2,-0.2,1.2]); return; pause(0.1)

  fU = fftn(U); 
  fU = dealias(fU);

  t = 0;
  n = 0;

%-- Main loop - time advancement --

  while t<tmax

    % plot and some diagnostics

    if (rem (n, dn) == 0)

      printf("n=%d,   t=%f\n", n, t);

      subplot(2,2,1);
      plot(x, U);  axis([-L/2,L/2,-0.2,1.2]);

      subplot(2,2,2);  axis([-N/2,N/2,1.e-3,10]);
      semilogy(wn(N), abs(fU), '*'); axis([-N/2,N/2,1.e-4, 1.e3]);

      subplot(2,2,3);
      plot(x, dU(U,L));  axis([-L/2,L/2,-0.1,0.1]);

      subplot(2,2,4);
      plot(x, ddU(U,L));  axis([-L/2,L/2,-0.02,0.02]);

      %return;
      pause(1.0);

    end

    %------------------------

    %fU = advance_splitstep(fU, L, dt);      U  =  ifftn(fU);
    fU = advance_euler_k(fU, L, dt);      U  =  ifftn(fU);
    %U = advance_euler_r(U, L, dt);
 
    t = t + dt;
    n = n + 1;

    %------------------------

  end


end

%-------------------------------------------------
function [U,x] = init_gauss(L, N, h0, r0)

  dx = L/N;

  x  = (-N/2:1:N/2-1)*dx + dx/2;
  U  = h0*exp( -(x/r0).^2 );

end


%-------------------------------------------------
function [U,x] = init_tanh(L, N, h0, r0)

  dx = L/N;

  x  = (-N/2:1:N/2-1)*dx + dx/2;

  n0 = 8*h0/dx;

  n  = 0:1:N-1;
  f1 = N/4;
  f2 = N*3/4;
  k  = 2.0/n0;

  f1 =  0.5*tanh((n-f1)*k);
  f2 = -0.5*tanh((n-f2)*k);

  U  = (f1 + f2)*h0;

end



%-------------------------------------------------
function fU = advance_splitstep(fUold,L,dt)
 
    N = length(fUold);
 
    DD = laplacian(N,L);
    Q = (1 + 0.25*dt*DD)./(1 - 0.25*dt*DD); 

    fU = fUold;

    % advance linear part half-step
    fU = Q .* fU; 
 
    % advance non-linear part whole step
    U  =  ifftn(fU);
    R  =  0.25*U.*(1-U);
    R  =  -R.*sign(dU(U,L));
    U  =  U + R*dt;
    fU =  fftn(U);
    fU =  dealias(fU);
 
    % advance linear part half-step
    fU = Q .* fU;
 
end

%-------------------------------------------------

function fU = advance_euler_k(fUold,L,dt)
 
    N = length(fUold);
 
    DD = laplacian(N,L);
    fU = fUold;

    % reaction
    U  =  ifftn(fU);
    R  =  0.25*U.*(1-U);
    R  =  -R.*sign(dU(U,L));
    fR =  fftn(R);
    fR =  dealias(fR);
 
    %fR = 0;

    % advancing solution in Fourier space
    fU = fU + (DD.*fU + fR)*dt;

end

%-------------------------------------------------

function Unew = advance_euler_r(U,L,dt)

    % reaction
    R  =  0.25*U.*(1-U);
    R  =  -R.*sign(dU(U,L));

    %R = 0;

    % diffusion

    D = ddU(U,L);

    % advancing solution in real space

    Unew  =  U + (D + R)*dt;

end


%-------------------------------------------------

function DD = laplacian(N,L)
% create a matrix of wavenumbers squared 
% normalized to the domain size

  k = wn(N);

  DD = -(k*2*pi/L).^2;

end


%-------------------------------------------------

function Ek=spectrum(U)

  N = length(U);

  fU = fftn(U) * 2*pi/N;

  Ek = abs(fU);

end

%-------------------------------------------------

function F=ddU(U,L)

  N = length(U);
  k = wn(N);

  DD = -(k*2*pi/L).^2;

  fU = fftn(U);
  fU = fU.*DD; 

  F  = ifftn(fU); 

end

%-------------------------------------------------

function F=dU(U,L)

  N = length(U);
  k = wn(N);

  D = i*k*2*pi/L;

  fU = fftn(U);
  fU = fU.*D; 

  F  = ifftn(fU);

end

%-------------------------------------------------

function k=wn(N)
% create the vector of k 
% k = [0,1,2...,N/2, -N/2, ..., -2, -1]

  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);

end
%-------------------------------------------------
%-------------------------------------------------

function fU=dealias(fUin)
% zero out high harmonics

  fU = fUin;

  N = length(fU);

  %fU(N/4:3*N/4+1) = 0;   % COMMENT / UNCOMMENT THIS LINE

end
%-------------------------------------------------
