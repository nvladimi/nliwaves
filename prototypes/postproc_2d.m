function postproc_2d(U, x, t, fbase)
%
% diagnostics of individual collapse 
%
% input:  U - 2D array of Psi
%         x - 1D array of coordiantes, x=y
%         t - time
%
% output: screen and files

  persistent fn;        % file number
  persistent H0;        % initial energy
  persistent NP0;       % initial number of particles


  %-- absolute value of Psi

  u = abs(U);
  dx = x(2)-x(1);

  %-- output file names -- 

  f_nph=strcat(fbase,'.nph');   % numberof particles and hamiltonian
  f_max=strcat(fbase,'.max');   % location, maximum, and width of collapse
  f_rad=strcat(fbase,'.rad');   % radial profile, u=u(r)
  f_fft=strcat(fbase,'.fft');   % energy spectra, Ek

  
  %-- initialize files and reference quantities --

  if t==0

    fn = 0;

    % compute intitial energy and number of particles

    [H1,H2] = hamiltonian(U,dx);
    H0 = H1-H2;
    NP0 = numparticles(U,dx);

    % write file headers

    outfile = fopen (f_nph, 'wt');
    fprintf(outfile, ...
            '#1.t   2.NP   3.H1   4.H2   5.errNP   6.errH \n\n'); 
    fclose(outfile);
  

    outfile = fopen (f_max, 'wt');
    fprintf(outfile, ...
            '#1.t  2.Psi_max  3.x_max  4.y_max  5.R_curv\n\n'); 
    fclose(outfile);

  end


  %-- Number of particles and hamiltonian --

  [H1,H2] = hamiltonian(U,dx);
  NP     = numparticles(U,dx);

  outfile = fopen (f_nph, 'at');
  fprintf(outfile, '%6.2f  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n', ...
           t, NP, H1, H2, (NP-NP0)/NP0, (H1-H2-H0)/H0); 
  fclose(outfile);


  %-- Location, maximum, and width of collapse  --

  [u0, x0, y0, c0] = collapse_max(u, x);

  outfile = fopen (f_max, 'at');
  fprintf(outfile, '%6.2f  %12.4e  %8.4f  %8.4f  %12.4e\n', ...
           t, u0, x0, y0, c0);
  fclose(outfile);


  %-- Profile and spectrum --

  x0=0; y0=0;
    
  fname = strcat(f_rad, num2str(fn,'.%04d'));
  collapse_profile(u, x, x0, y0, u0, fname, t);

  fname = strcat(f_fft, num2str(fn,'.%04d'));
  collapse_spectrum(U, fname, t);


  %-- Prepare for the next output --
 
  fn = fn + 1;

end

%-------------------------------------------------

function NP = numparticles(U,dx)
% compute number of particles

  NP = sum(sum(abs(U).^2))*dx^2;

end
%-------------------------------------------------

function [H1,H2] = hamiltonian(U,dx)
% compute two terms of hamiltonian, H = H1 - H2

  N  = length(U);

  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);

  [k1, k2] = meshgrid (k, k);

  kk = k1.^2 + k2.^2;

  fU = fftn(U);

  H1 = (2*pi/N^2)^2 * sum(sum(kk.*abs(fU).^2));

  H2 = 0.5*dx^2 * sum(sum(abs(U).^4));

end

%-------------------------------------------------

function collapse_profile(u,x,x0,y0,u0,fname,t)

  N = length(u);
  y = x;

  [xx, yy] = meshgrid (x-x0, y-y0);

  r  = sqrt(xx.^2 + yy.^2);
  r  = reshape(r,N*N,1);
  ur = reshape(u,N*N,1);

  [r,ind] = sort(r);
  ur = ur(ind);

  uavg  = sum(sum(u))/(N*N);

  % write profile to the file 

  outfile = fopen (fname, 'w');
  fprintf(outfile, '#1.r  2.abs(Psi) \n' );
  fprintf(outfile, '#time = %12.8f\n\n', t);
  %fprintf(outfile, '%12e  %12e \n', 0, u0);
  for i=1:N
     fprintf(outfile, '%12e  %12e \n', r(i), ur(i));
     if ur(i)<uavg 
        break;
     end
  end
  fclose(outfile);

end
%-------------------------------------------------

function collapse_spectrum(U,fname,t)

  N = length(U);

% create the vector of k 
% k = [0,1,2...,N/2, -N/2, ..., -2, -1]
   
  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);

  [k1, k2] = meshgrid(k, k);

  kk = round(sqrt(k1.^2 + k2.^2));

  fU = fftn(U) * 2*pi/N^2;

  fU = (abs(fU)).^2;
 
  Ek = zeros(N,1);

  for k1=1:N
    for k2=1:N

       k = kk(k1,k2)+1;
       if (k < N) 
         Ek(k) = Ek(k) + fU(k1,k2);
       end

    end
  end

  % write spectrum to the file

  outfile = fopen (fname, 'wt');
  fprintf(outfile, '#1.k  2.Ek \n\n' );
  fprintf(outfile, '#time = %12.8f\n\n', t);
  for i=1:N
     fprintf(outfile, '%12e  %12e \n', i-1, Ek(i));
  end
  fclose(outfile);

end

%-------------------------------------------------

function [u0,x0,y0,c0] = collapse_max(u,x)
% location, height, and width of individual collapse
%
% input:  u - 2D array of |psi|
%         x - 1D array of coordiantes, x=y
%
% output: [ psi_max, x_center, y_center, R_curvature ]

  y = x;      % grid is assumed to be square and uniform  

  m = 1;      % number of points around maximum used for interpolation
  n = 2*m+1;  % interpoalte over nxn square

  % Find rough (index) location of maximum.

  [uj, jj]  = max(u);
  [umax,i0] = max(uj);
  j0        = jj(i0);

  % For smooth location, interpolate Psi over nxn points near the maximum with
  % u(x,y) = c1 + c2*x + c3*y - c4*(x^2+y^2).
  % To find coefficients solve linear system Ac=F, where A=A(x,y), F=u(x,y).

  A = zeros(n^2,4);
  F = zeros(n^2,1);
 
  k = 0;

  for i=i0-m:i0+m
    for j=j0-m:j0+m
      k      = k+1;
      A(k,1) =  1;
      A(k,2) = x(i);
      A(k,3) = y(j);
      A(k,4) = -(x(i)^2 + y(j)^2);
      F(k)   =  u(j,i);
    end
  end

  c = A \ F;
 
  % coefficients in u(x,y) = u0 - c0*[ (x-x0)^2 + (y-y0)^2 ]
  
  c0 = c(4);
  x0 = c(2)/(2*c0);
  y0 = c(3)/(2*c0);
  u0 = c(1) + c0*(x0^2 + y0^2);

  % make c0 the radius of curvature at the apex

  c0 = 1./(2*c0);
  
end

%-------------------------------------------------
