function dse_fixed_point

%-- Input parameters --  

  outfile = 'grounsol_test.dat'   %'groundsol_nu05_rho-1_1024.dat'

  nu = 0.1;
  rho =  -1.0;
  alpha = 3/2;


  L     =   25.6;           % length of the box
  r0    =  sqrt(2);         % initial width of gaussian psi
  h0    =    1.0;           % initial height of gaussian psi
  x0    =    0.0;           % initial x-location of gaussian psi
  y0    =    0.0;           % initial y-location of gaussian psi      

  N     =   1024;           % discretization points 

  dx    = L/N;

%-- Initial data --
 
  [F, x]  = initPsi_2D(L, N, h0, r0, x0, y0);
  %imagesc (F); axis('equal', 'tics', 'off'); return;


%-- Grid and multipliers in k-space --

% k = [0,1,2...,N/2, -N/2, ..., -2, -1]
   
  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]) * (2*pi/L);

  [kx, ky] = meshgrid (k, k);

  Mf = 1 ./ (1 + kx.*kx + ky.*ky) ;
  Mv = kx.*kx  ./ (kx.*kx + nu*ky.*ky); 
  Mv(1,1) = 0;


%-- Fixed point iterations --

for s=1:100 

    fV = Mv .* fft2( F.*F );

    V = real(ifft2(fV));
    R = (F.*F - rho*V) .* F;

    fR = fft2(R);
    fF = fft2(F);
   
    SR = Mf .* fR .* conj(fF) ; 
    SL = fF .* conj(fF); 

    factor = ( real(sum(sum(SL))) / real(sum(sum(SR))) ) ^ alpha;

    fF = Mf .* fR * factor ;

    Fnew = real (ifft2(fF));

    Err = max(max(Fnew - F));
    F = Fnew;

    Fmax =  max(max(F));
    %min(min(F))
    disp([s, Fmax, Err, factor]);

    %imagesc(F); axis('equal', 'tics', 'off');

end
 
fV = Mv .* fft2( F.*F );
V  = real(ifft2(fV));
Vo =  V(N/2+1, N/2+1);
   
imagesc(V); axis('equal', 'tics', 'off');



%-- astigmatism --

  fF = fft2(F.*F);

  Fx = sum(sum(abs(ifft2(i*kx.*fF))));
  Fy = sum(sum(abs(ifft2(i*ky.*fF))));

  Fxx = real(ifft2(-kx.*kx.*fft2(F)));
  Fyy = real(ifft2(-ky.*ky.*fft2(F)));

  %max(max(F)), F(N/2+1, N/2+1)
  %imagesc(Fxx); ; axis('equal', 'tics', 'off'); max(max(abs(Fxx)))

  format long 

  h = abs(F(N/2+1, N/2+1));
  hxx = Fxx(N/2+1, N/2+1);
  hyy = Fyy(N/2+1, N/2+1);

  deriv_ratio = hxx / hyy;

  astigmatism = Fx/Fy;

  num_waves  = sum(sum(F.*F))*dx*dx;

  disp([nu, rho, num_waves, h, hxx, hyy, astigmatism, Vo])

  return

    %-- save binary F(x,y) --

  fid = fopen(outfile, 'wb');
  fwrite(fid, F, 'double');
  fclose(fid);

end


%-------------------------------------------------
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
