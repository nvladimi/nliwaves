function naive
% solves 1-D NLS equation 
%
%   i dPsi/dt + (1-i*a*eps)*del^2(Psi) + (1+i*eps)*|Psi|^2*Psi = i*b*eps*Psi 
%
% using naive central differences

% Input parameters  

  coefA =  1;               % coefficient "a" in equation 
  coefB =  0;               % coefficient "b" in equation 
  coefE =  0.01;            % coefficient "epsilon"

  L     = 12.8;             % length of the box
  r0    =  1.0;             % initial width of gaussian psi
  h0    =  3.6;             % initial height of gaussian psi
  tmax  =  0.213;           % maximum time

  N     =  8192;            % discretization points 
  fbase = 'r1h36';          % name base for output files
  nout1 =   256;            % save maximum every (nout1) iterations
  nout3 =   40960;          % save solution every (nout3) iterations 
  
  CFL   = 0.01;

%--- Derived parameters ---

  c1    = i + coefA*coefE;  % complex coefficients 
  c2    = i - coefE;
  c3    = coefB*coefE;

%--- Initial conditions ---

  [S,x] = initPsi(L, N, h0, r0);

  t     = 0;
  iter  = 0;
  postproc_1r(S,x,t,fbase,0,3);
  postproc_1r(S,x,t,fbase,0,1);

%--- Evolve ---

  dx    = L/N;                                     % grid spacing 
  d0=55/24; d1=-59/24; d2=37/24; d3=-9/24;         % adam bashford 4th coefficients

  while t<tmax

    dt = CFL*dx*dx;

    if t~=0                                        % rotate RHSs from previous steps
       s3=s2; s2=s1; s1=s0;
       dt3=dt2; dt2=dt1; dt1=dt0; dt0=dt;
    end

    DD  = d2r(S, x, dx);                           % laplacian in radial coords
    s0 = c1*DD + (c2*abs(S).^2 + c3).*S;           % RHS of the equation

    if t==0                                        % no saved RHSs at t=0, 
       s1=s0; s2=s0; s3=s0;                        % use current
       dt3=dt; dt2=dt; dt1=dt; dt0=dt;
    end

    S = S + s0*d0*dt0 + s1*d1*dt1 + s2*d2*dt2 + s3*d3*dt3;  % evolve solution
    t = t + dt0;
    iter = iter + 1;

    %-- output --

    if rem(iter, nout1) == 0
      %plot(x, abs(S), 'xr');  axis([-L/2,L/2,0,16]); pause(0.1);
      %grid on; %return

      fn = iter/nout1;
      postproc_1r(S,x,t,fbase,fn,1);
    end

    if rem(iter, nout3) == 0
      fn = iter/nout3;
      postproc_1r(S,x,t,fbase,fn,3);
    end 

  end


end
%-------------------------------------------------
function [U,x] = initPsi(L, N, h0, r0)

  dx = L/N;

  x  = (-N/2:1:N/2-1)*dx + dx/2;
  U  = h0*exp( -(x/r0).^2 ) + 0*i;

end

%--------------------------------------------------

function f = d2x(u,dx)
% compute second derivative with central differences, 
% periodic boundary conditions

  k = 1/(dx*dx*12);
  f = (-circshift(u,[0,2]) + 16*circshift(u,[0,1]) - 30*u + ...
       16*circshift(u,[0,-1]) - circshift(u,[0,-2]) )*k;

end

%--------------------------------------------------

function f = d2r(u,x,dx)
% compute laplacian in radial coords with central differences, 
% periodic boundary conditions (kinda silly)

  k2 = 1/(dx*dx*12); 
  k1 = 1/(dx*12);

  d2 = (-circshift(u,[0,2]) + 16*circshift(u,[0,1]) - 30*u + ...
       16*circshift(u,[0,-1]) - circshift(u,[0,-2]) )*k2;

  d1 = (circshift(u,[0,2]) - 8*circshift(u,[0,1]) + ...
       8*circshift(u,[0,-1]) - circshift(u,[0,-2]) )*k1;

  d1 = d1./x;

  f = d1 + d2;

end


%-------------------------------------------------

