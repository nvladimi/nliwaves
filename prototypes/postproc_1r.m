function postproc_1r(U, x, t, fbase, fn, task)
%
% diagnostics of individual collapse 
%
% input:  U     - 1D array of Psi(r)
%         x     - 1D array of coordiantes, x=r
%         t     - time
%         fbase - file name base
%         fn    - output number
%         task  - "what to do" flag
%
% output: text files

  switch (task)
     case 1
       collapse_max(U, x, t, fbase, fn);
     case 2
       hamiltonian(U, x, t, fbase, fn);
     case 3
       collapse_profile(U, x, t, fbase, fn);
     otherwise % do all tasks
       collapse_max(U, x, t, fbase, fn);
       hamiltonian(U, x, t, fbase, fn);
       collapse_profile(U, x, t, fbase, fn);
   end
 
end

%-------------------------------------------------

function hamiltonian(U, x, t, fbase, fn)
% compute and  write to file number of particles and hamiltonian

  persistent H0;        % initial energy
  persistent NP0;       % initial number of particles

  fname=strcat(fbase,'.nph');
  dx = x(2)-x(1);

  %-- compute number of particles --

  NP = sum(sum(abs(U).^2).*x)*dx*2*pi;

  %-- compute two terms of hamiltonian, H = H1 - H2 --

  dU = (circshift(U,[0,2]) - 8*circshift(U,[0,1]) + ...
       8*circshift(U,[0,-1]) - circshift(U,[0,-2]) )/(dx*12);

  dU(1)     = (-3*U(1)+4*U(2)-U(3))/(2*dx); 
  dU(2)     = (U(3)-U(1))/(2*dx); 
  dU(end-1) = (U(end)-U(end-2))/(2*dx); 
  dU(end)   = (3*U(end)-4*U(end-1)+U(end-2))/(2*dx); 

  H1 =     sum(sum(abs(dU).^2 .* x))*2*pi*dx;

  H2 = 0.5*sum(sum(abs(U).^4 .* x))*2*pi*dx;

  %-- write text data  --


  if fn==0

    H0 = H1-H2;
    NP0 = NP;

    fid = fopen(fname, 'wt');
    fprintf(fid, '1.t   2.NP   3.H1   4.H2   5.errNP   6.errH \n\n'); 
    fclose(fid);

  end


  fid = fopen (fname, 'a');
  fprintf(fid, '%12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n', ...
           t, NP, H1, H2, (NP-NP0)/NP0, (H1-H2-H0)/H0); 
  fclose(fid);

end

%-------------------------------------------------

  function collapse_profile(U, x, t, fbase, fn)
% profile of collapse, U(r) 

  N = length(U);

  fbase=[fbase,'.rad'];
fname = strcat(fbase, num2str(fn,'.%04d'));

  fid = fopen (fname, 'wt');
  fprintf(fid, '#1.r  2.Re(Psi)  3.Im(Psi)  4.abs(Psi) \n' );
  fprintf(fid, '#time = %12.8f\n\n', t);

  for i=N/2+1:N
     fprintf(fid, '%12.4e  %12.4e  %12.4e  %12.4e\n', ...
	     x(i), real(U(i)), imag(U(i)), abs(U(i)) );
  end
  fclose(fid);

end

%-------------------------------------------------

function u0 = collapse_max(U, x, t, fbase, fn)
% height and width of individual collapse

  fname=strcat(fbase,'.max');

  %-- Find index location of maximum --

  N  = length(U);
  u  = U(N/2-1:N/2+2);   % select four center points
  dx = x(2)-x(1); 

  %-- Find derivatives --

  d0 = ( - u(1) +  9*u(2) +  9*u(3) - u(4) ) / 16;
  d1 = (   u(1) - 27*u(2) + 27*u(3) - u(4) ) /(24*dx);
  d2 = (   u(1) -    u(2) -    u(3) + u(4) ) /(2*dx^2);

  h   = abs(d0);
  u   = real(d0);
  v   = imag(d0);
  p   = atan2(v,u);
  ddu = real(d2);
  ddv = imag(d2);
  ddh = (u*ddu + v*ddv)/h;
  ddp = (u*ddv - v*ddu)/(h*h);
   
  %-- write data to file --

  if fn==0
    fid = fopen (fname, 'wt');
    fprintf(fid, '#1.t  2.h  3.ddh  4.ddp  5.empty  6.empty  7.empty  8.empty\n');
    fprintf(fid, '#9.empty   10.empty  11.phase  12.u  13.v  14.ddu  15.ddv\n');
    fprintf(fid, 'NaN NaN NaN NaN NaN  NaN NaN NaN NaN NaN  NaN NaN NaN NaN NaN\n\n');
    fclose(fid);
  end

  fmt='%20.12e %16.8e %16.8e %16.8e  0 0 0 0 0 0  %16.8e %16.8e %16.8e %16.8e %16.8e\n';

  fid = fopen (fname, 'at');
  fprintf(fid, fmt,  t, h, ddh, ddp,   p, u, v, ddu, ddv);
  fclose(fid);
 
end

%-------------------------------------------------
