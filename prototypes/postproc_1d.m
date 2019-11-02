function postproc_1d(U, x, t, fbase)
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

    outfile = fopen(f_nph, 'wt');
    fprintf(outfile, '1.t   2.NP   3.H1   4.H2   5.errNP   6.errH \n\n'); 
    fclose(outfile);
  

    outfile = fopen (f_max, 'wt');
    fprintf(outfile, '#1.t  2.Psi_max  3.x_max  4.y_max  5.R_curv\n\n'); 
    fclose(outfile);

  end


  %-- Number of particles and hamiltonian --

  [H1,H2] = hamiltonian(U,dx);
  NP     = numparticles(U,dx);

  outfile = fopen (f_nph, 'a');
  fprintf(outfile, '%6.2f  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n', ...
           t, NP, H1, H2, (NP-NP0)/NP0, (H1-H2-H0)/H0); 
  fclose(outfile);


  %-- Maximum and width of collapse  --

  [u0, c0] = collapse_max(u, x);

  outfile = fopen (f_max, 'at');
  fprintf(outfile, '%6.2f  %12.4e  %12.4e\n', t, u0, c0);
  fclose(outfile);

  %-- Profile and spectrum --
    
  fname = strcat(f_rad, num2str(fn,'.%04d'));
  collapse_profile(u, x, u0, fname);

  %fname = strcat(f_fft, num2str(fn,'.%04d'));
  %collapse_spectrum(U, fname);


  %-- Prepare for the next output --
 
  fn = fn + 1;

end

%-------------------------------------------------

function NP = numparticles(U,dx)
% compute number of particles

  NP = sum(sum(abs(U).^2))*dx;

end
%-------------------------------------------------

function [H1,H2] = hamiltonian(U,dx)
% compute two terms of hamiltonian, H = H1 - H2

  N  = length(U);
  L  = N*dx;

  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);

  fU = fftn(U);

  kk = k.^2;


  H1 = (L/N^2)*(2*pi/L)^2 * sum(sum(kk.*abs(fU).^2));

  H2 = 0.5*dx * sum(sum(abs(U).^4));

end

%-------------------------------------------------

function collapse_profile(u,x,u0,fname)

  N = length(u);

  % write profile to the file 

  outfile = fopen (fname, 'wt');
  fprintf(outfile, '#1.r  2.abs(Psi) \n\n' );
  %fprintf(outfile, '%12e  %12e \n', 0, u0);
  for i=N/2+1:N
     fprintf(outfile, '%12e  %12e \n', x(i), u(i));
  end
  fclose(outfile);

end

%-------------------------------------------------

function collapse_spectrum(U,fname)

  N = length(U);

% create the vector of k 
% k = [0,1,2...,N/2, -N/2, ..., -2, -1]
   
  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);

  fU = fftn(U) * 2*pi/N;

  Ek = (abs(fU)).^2;
 
  % write spectrum to the file

  outfile = fopen (fname, 'wt');
  fprintf(outfile, '#1.k  2.Ek \n\n' );
  for i=1:N/2
     fprintf(outfile, '%12e  %12e \n', k(i), Ek(i));
  end
  fclose(outfile);

end

%-------------------------------------------------

function [u0,c0] = collapse_max(u,x)
% height and width of individual collapse
% maximum is assumed to be in the center
%
% input:  u - 2D array of |psi|
%         x - 1D array of coordiantes, x=y
%
% output: [ psi_max, R_curvature ]

  % Find index location of maximum.

  N = length(u);
  n1 = N/2+1;
  n2 = n1+1;

  % Find parabola passing throgh two points

  x1 = x(n1);
  x2 = x(n2);

  u1 = u(n1);
  u2 = u(n2);

  c0 = (u1-u2)/(x2^2 - x1^2);
  u0 = u1+c0*x1^2;


  % make c0 the radius of curvature at the apex

  c0 = 1./(2*c0);
  
end

%-------------------------------------------------
