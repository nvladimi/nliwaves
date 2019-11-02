function gain = linSE(fbase, N, L, t_all)
% solves NLS equation 
%
%   i dPsi/dt + del^2(Psi) = 0 
%
% using 2nd order split-step method

%-- Initializaion --

  iunit = i;
 
  U = read_psi([fbase, '.psi.0000'], N);

  max_psi_0 = max(max(abs(U)));

  h_out =[];

%-- precompute exponents for linear step --

  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);
  [kx,ky] = meshgrid(k,k);

%-- time advancement --

  DD =  -(kx.^2 + ky.^2)*(2*pi/L)^2;

  fU0 = fftn(U);            % FT to k-space

  for t = t_all

    expLin = exp(iunit*DD*t);

    fU = fU0 .* expLin;     % time advancement
    U  = ifftn(fU);         % FT to r-space

    max_psi = max(max(abs(U)));

    h_out = [h_out; t, max_psi]; 

  end

  figure(1); imagesc(abs(U), [0,3]);     axis("equal", "off");
  figure(2); imagesc(arg(U), [-pi,pi]);  axis("equal", "off");

  max_psi = max(max(abs(U)));

  gain = max_psi / max_psi_0;

  save([fbase, '.lin.dat'], 'h_out');

end


%-------------------------------------------------

