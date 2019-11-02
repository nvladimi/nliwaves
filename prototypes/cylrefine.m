function cylrefine

% solves 1-D NLS equation 
%
%   i dPsi/dt + (1-i*a*eps)*del^2(Psi) + (1+i*eps)*|Psi|^2*Psi = i*b*eps*Psi 
%
% using naive central differences

% -- Input parameters --------------------------------------  

  coef(1) =  1;               % coefficient "a" in equation 
  coef(2) =  0;               % coefficient "b" in equation 
  coef(3) =  0.01;            % coefficient "epsilon"

  L       =   12.8;           % length of the box
  r0      =    1.0;           % initial width of gaussian psi
  h0      =    3.4;           % initial height of gaussian psi
  tmax    =    0.263;         % maximum time
  nout1   =    8;             % save maximum every (nout1) iterations
  nout3   =    160;           % save solution every (nout3) iterations 

  N0      =  512;             % number of coarse cells
  N1      =  256;             % initial number of refined cells
  RF      =     2;            % initial refinement factor
  fbase0  = 'r1h34a4';        % output files for the coarse data
  fbase1  = 'r1h34b4';        % output files for the fine data

  CFL     =    0.01;


%-- Initial conditions --------------------------------------

  [U0,x0] = initPsi(L,       N0,    h0, r0);
  [U1,x1] = initPsi(L*N1/N0, N1*RF, h0, r0);

  t       = 0;
  iter    = 0;

  postproc_1r(U0,x0,t,fbase0, 0, 1);
  postproc_1r(U1,x1,t,fbase1, 0, 1);
  postproc_1r(U0,x0,t,fbase0, 0, 3);
  postproc_1r(U1,x1,t,fbase1, 0, 3);


%-- Setup mesh ----------------------------------------------

  dx0  = L/N0;                % coarse grid spacing
  dx1  = dx0/RF;              % fine grid spacing
  dt0  = CFL*dx0*dx0;         % coarse timestep

  ind1 = (RF/2):RF:(N1*RF);   % fine points left  of coarse points
  ind2 = ind1 + 1;            % fine points right of coarse points
  i1   = N0/2 - N1/2 + 1;     % most left  coarse point on the fine mesh
  i2   = N0/2 + N1/2;         % most right coarse point on the fine mesh

  Uref = 2*max(abs(U0));      % refine when Umax > Uref

%-- Evolve solution -----------------------------------------

  while t<tmax

    %-- refine mesh if needed --

    if ( max(abs(U1)) > Uref && N1>4)

       RF   = RF*2;
       N1   = N1/2;
       Uref = Uref*2;

       ind1 = (RF/2):RF:(N1*RF);
       ind2 = ind1+1;
       i1   = N0/2 - N1/2 + 1;
       i2   = N0/2 + N1/2;

       dx1  = dx0/RF;

       x1   = refine(x1);
       U1   = refine(U1);

       is_refine = 1;

       nout1 = nout1/2;

    else 

       is_refine = 0;

    end

    %-- fine timestep --
    dt1  = dt0/RF^2;

    %-- evolve coarse solution --
    U0 = U0 + dUcoarse(U0,x0,dt0,coef);

    %-- evolve fine solution --
    for k=1:RF^2
        U1 = U1 + dUfine(U1,x1,dt1,coef,is_refine);
    end

    %-- adjust fine data at the left merge --
    d0 = dEdgeC( U0(i1-2:i1+1), dx0 );          % coarse data and derivative
    d1 = dEdgeL( U1(1:3), dx1 );                % fine   data and derivative
    dU = adjustL( d0-d1, dx0, RF );             % interpolate difference
    U1(1:RF) = U1(1:RF) + dU;                   % adjust data

    %-- adjust fine data at the right merge --
    d0 = dEdgeC( U0(i2-1:i2+2), dx0 );          % coarse data and derivative
    d1 = dEdgeR( U1(end-2:end), dx1 );          % fine   data and derivative
    dU = adjustR( d0-d1, dx0, RF );             % interpolate difference
    U1(end-RF+1:end) =  U1(end-RF+1:end) + dU;  % adjust data

    %-- update coarse solution with fine data --
    U0(i1:i2) = 0.5*( U1(ind1) + U1(ind2) );

    t = t + dt0;
    iter = iter + 1;

    %-- output --
    if rem(iter, nout1) == 0
      %plot(x0, abs(U0), 'xr', x1, abs(U1),'-b');  axis([-L/2,L/2,0,16]); pause(0.1);
      %grid on; %return

      fn = iter/nout1;
      postproc_1r(U0,x0,t,fbase0,fn,1);
      postproc_1r(U1,x1,t,fbase1,fn,1);
    end

    if rem(iter, nout3) == 0
	fn = iter/nout3;
        postproc_1r(U0,x0,t,fbase0,fn,3);
        postproc_1r(U1,x1,t,fbase1,fn,3);
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

function dU = adjustL(d, L, n)

    % coefficients in cubic polynomial

    L = -L;

    a = (3*d(1) - d(2)*L) / L^2;
    b = (d(2)*L - 2*d(1)) / L^3;

    dx = L/n;
    x  = dx/2:dx:L;

    dU = a*(L-x).^2 + b*(L-x).^3;

end

%--------------------------------------------------

function dU = adjustR(d, L, n)

    % coefficients in cubic polynomial

    a = (3*d(1) - d(2)*L) / L^2;
    b = (d(2)*L - 2*d(1)) / L^3;

    dx = L/n;
    x  = dx/2:dx:L;

    dU = a*x.^2 + b*x.^3;

end

%--------------------------------------------------

function d = dEdgeC(U,dx)

    d(1) = ( - U(1) +  9*U(2) +  9*U(3) - U(4) ) / 16;
    d(2) = (   U(1) - 27*U(2) + 27*U(3) - U(4) ) /(24*dx);
    d(3) = (   U(1) -    U(2) -    U(3) + U(4) ) /(2*dx^2);

end

function d = dEdgeL(U,dx)

    d(1) = ( 15*U(1) - 10*U(2) + 3*U(3)) / 8;
    d(2) = ( -2*U(1) +  3*U(2) -   U(3)) / dx;
    d(3) = (    U(1) -  2*U(2) +   U(3)) /(dx^2);

end

function d = dEdgeR(U,dx)

    d(1) = ( 15*U(3) - 10*U(2) + 3*U(1)) / 8;
    d(2) = (  2*U(3) -  3*U(2) +   U(1)) / dx;
    d(3) = (    U(3) -  2*U(2) +   U(1)) /(dx^2);

end

%--------------------------------------------------

function dU = dUcoarse(U,x,dt,coef)

  persistent u0;
  persistent u1;
  persistent u2;
  persistent u3;
  persistent dt0;
  persistent dt1;
  persistent dt2;
  persistent dt3;
  
  c1 = i + coef(1)*coef(3);                            % complex coefficients 
  c2 = i - coef(3);
  c3 = coef(2)*coef(3);

  dx = x(2)-x(1);

  d0=55/24; d1=-59/24; d2=37/24; d3=-9/24;             % adams-bashforth coefficients

  if ~isempty(u0)                                      % rotate RHSs
    u3=u2; u2=u1; u1=u0;                               % from previos step
    dt3=dt2; dt2=dt1; dt1=dt0; dt0=dt;
  end

  DD  = d2r(U, x, dx);                                 % laplacian in radial coords

  u0 = c1*DD + (c2*abs(U).^2 + c3).*U;                 % RHS of the equation

  if isempty(u1)                                       % no saved RHSs at t=0, 
     u1=u0; u2=u0; u3=u0;                              % use current
     dt3=dt; dt2=dt; dt1=dt; dt0=dt;
  end

  dU = u0*d0*dt0 + u1*d1*dt1 + u2*d2*dt2 + u3*d3*dt3;  % change in solution

end

%--------------------------------------------------

function dU = dUfine(U,x,dt,coef,is_refine)

  persistent u0;
  persistent u1;
  persistent u2;
  persistent u3;
  persistent dt0;
  persistent dt1;
  persistent dt2;
  persistent dt3;

  d0=55/24; d1=-59/24; d2=37/24; d3=-9/24;             % adams-bashforth coefficients

  c1 = i + coef(1)*coef(3);                            % complex coefficients 
  c2 = i - coef(3);                                    %    in the NLS equation
  c3 = coef(2)*coef(3);

  if ~isempty(u0);                                     % rotate RHSs
    u3=u2; u2=u1; u1=u0;                               %    from previos steps
    dt3=dt2; dt2=dt1; dt1=dt0; dt0=dt;
  end

  if is_refine
     u1  = refine(u1);
     u2  = refine(u2);
     u3  = refine(u3);
  end

  dx    = x(2)-x(1);
  DD    = d2r(U, x, dx);                              % laplacian in radial coords

  ddU1  = (2*U(1)-5*U(2)+4*U(3)-U(4))/dx^2;           % one-sided derivatives
  ddU2  = (U(1)-2*U(2)+U(3))/dx^2;                    %    at the edges
  ddU3  = (U(end)-2*U(end-1)+U(end-2))/dx^2;
  ddU4  = (2*U(end)-5*U(end-1)+4*U(end-2)-U(end-3))/dx^2; 

  dU1   = (-3*U(1)+4*U(2)-U(3))/(2*dx); 
  dU2   = (U(3)-U(1))/(2*dx); 
  dU3   = (U(end)-U(end-2))/(2*dx); 
  dU4   = (3*U(end)-4*U(end-1)+U(end-2))/(2*dx); 

  DD(1)     = ddU1 + dU1/x(1);
  DD(2)     = ddU2 + dU2/x(2);
  DD(end-1) = ddU3 + dU3/x(end-1);
  DD(end)   = ddU4 + dU4/x(end);

  u0 = c1*DD + (c2*abs(U).^2 + c3).*U;                 % RHS of the equation

  if isempty(u1);                                      % no saved RHSs at t=0, 
     u1=u0; u2=u0; u3=u0;                              % use current
     dt3=dt; dt2=dt; dt1=dt; dt0=dt;
  end

  dU = u0*d0*dt0 + u1*d1*dt1 + u2*d2*dt2 + u3*d3*dt3;  % change in solution

end

%--------------------------------------------------

function Uout = refine(U)
% very clumsy... need to be rewritten

  N=length(U);

  Uc = ( -circshift(U,[0,1]) + 9*U + ...
       9*circshift(U,[0,-1]) - circshift(U,[0,-2]) )/16;

  U(1:2:2*N)=U;
  U(2:2:2*N)=Uc;

  Uc = ( -circshift(U,[0,1]) + 9*U + ...
       9*circshift(U,[0,-1]) - circshift(U,[0,-2]) )/16;

  Uout = Uc(N/2:3*N/2-1);

end

%--------------------------------------------------

function f = d2r(u,x,dx)
% compute laplacian in radial coords with central differences, 
% periodic boundary conditions

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

