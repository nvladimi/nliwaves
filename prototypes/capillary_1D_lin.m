function lin_profile
%
% Sinusoidal wave with binary output for 2D code.
% 

fbase = 'lin_A01_N256'

N = 256             % number of output points
L = 4*pi            % domain size
k = 2*pi/L * 2      % two wavelengths
A = 0.1
  
iunit = i;

%----------------------------------------------------

  x = (0:N-1)*L/N;

  psi = 1 + A*cos(k*x) - iunit * A*sin(k*x);

  plot(x,real(psi), x,imag(psi));  set(gca, "fontsize", 20);

  psi = repmat(psi, N, 1);

  p  = zeros(2,N,N);
  p(1,:,:) = real(psi);
  p(2,:,:) = imag(psi);

  p = reshape(p, N*N*2, 1);

  fid = fopen([fbase, '.psi.0000'], 'wb');
  fwrite(fid, p, 'double');
  fclose(fid);


end

%---------------------
