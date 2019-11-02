function test_fft

%-- Input parameters --  


  L     =   12.8;           % length of the box
  r0    =    1.0;           % initial width of gaussian psi
  h0    =    4.0;           % initial height of gaussian psi
  x0    =    0.0;           % initial x-location of gaussian psi
  y0    =    0.0;           % initial y-location of gaussian psi      

  N     =    256;           % discretization points 

  dx    = L/N;

%-- Initial data --
 
  [U, x]  = initPsi_2D(L, N, h0, r0, x0, y0);
  imagesc (real(U));


%-- dU/dx --

  V = Ux(U, L);
  imagesc (real(V));


%-- d2U/dx2

  V = Uxx(U, L);
  imagesc (real(V));


%-- Laplacian

  V = ddU(U, L);
  imagesc (real(V));

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

function DD = laplacian(N,L)
% create a matrix of wavenumbers squared 
% normalized to the domain size
% create the vector of k 
% k = [0,1,2...,N/2, -N/2, ..., -2, -1]
   
  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);


  [k1, k2] = meshgrid (k, k);
  DD = -(k1.^2 + k2.^2)*(2*pi/L)^2;

end

%-------------------------------------------------

function F=ddU(U,L)

  N = length(U);

  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);


  [k1, k2] = meshgrid (k, k);
  DD = -(k1.^2 + k2.^2)*(2*pi/L)^2;

  fU = fft2(U);
  fU = fU.*DD; 

  F  = ifft2(fU); 

end

%-------------------------------------------------

function F=Uxx(U,L)

  N = length(U);

  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);

  [k1, k2] = meshgrid (k, k);
  DD = -(k1.^2)*(2*pi/L)^2;

  fU = fft2(U);
  fU = fU.*DD; 

  F  = ifft2(fU); 

end


%-------------------------------------------------

function F=Ux(U,L)

  N = length(U);

  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);

  [k1, k2] = meshgrid (k, k);
  DD = i*k1*(2*pi/L)^2;

  fU = fft2(U);
  fU = fU.*DD; 

  F  = ifft2(fU); 

end


%-------------------------------------------------
