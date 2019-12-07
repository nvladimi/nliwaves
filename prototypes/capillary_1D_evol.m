function capillary_1D_evol
%
%  Evolution of 1D capillary wave in velocity formulation (unfinished).
%  Need good input profile for elevation with clean 2nd derivative.
%  Currently only input profile and derivatives are plotted. 
%



%-- Input parameters --  

  L = 1.2656000000e+01 
  N = 256;
  dt = 0.01;

  fname = 'a1_c1p01_b-1.txt';

%-------------------

  data = load(fname);

  h = transpose(data(:,2));
  u = transpose(data(:,3));


  dx    = L/N;

%  x  = (-N/2:1:N/2-1)*dx + dx/2;
  x  = (0:N-1)*dx;
  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k1 = circshift(k-k0,[0,-k0]);

  DD = -(k1.^2)*(2*pi/L);
  D  = i*k1*(2*pi/L);


%-- Test derivatives --
 
 
  U = h;
figure(1); plot(x, U, '-o'); set(gca, "fontsize", 20);

  V = real(ifft(fft(U).*D)); 
  figure(2); plot(x, V); set(gca, "fontsize", 20);

  V  = real(ifft(fft(U).*DD)); 
figure(3); plot(x, V, '-o'); set(gca, "fontsize", 20);

return

%-- Evolution --

  t=0;

  for n=1:1

    hxx = real(ifft(fft(h).*DD));

    f1 = - h.*u;
    f2 = hxx - 0.5*u.*u;

    dh = real(ifft(fft(f1).*D));
    du = real(ifft(fft(f2).*D)); 

    h = h + dh*dt;
    u = u + du*dt;
    t = t + dt;

  end

  
  figure(2); plot(x, h); set(gca, "fontsize", 20);
  figure(3); plot(x, u); set(gca, "fontsize", 20);


  
end
