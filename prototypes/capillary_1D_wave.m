function capillary_1D_wave
%
%  An attempt to integrate ODE to obtain a traveling wave profile (unfinished).
%  Currently, the solution have issues with second derivative at x=0 and x=lambda/2.
%  Need more accurate stitching of half-wave-lenghts.
%  Parameters "b" and "c" must be selected to provide periodic potential (<u> = 0). 
%


global a; global b; global c;

fbase = 'a1_c1p01_b-1'

a = 1      % wave speed
c = 1.01   % controls amplitude
b = -1      % controls inflection point

N = 256      % number of output points

tmax = 4;
dt = 0.0001;

iunit = i;

%----------------------------------------------------

%-- profile height limits

z1 = c - sqrt(c*c - b*b/a/a);
z2 = c + sqrt(c*c - b*b/a/a);

%-- intergrate half-wavelenght forward and backward

   t2 = 0:dt:tmax;
   t1 = -t2;

   z0 = c;
   p0 = a * (1 - 1/z0);
   
   f0 = [z0, p0];

   f1 = lsode('fdot', f0, t1);
   f2 = lsode('fdot', f0, t2);

%-- trim to min and max and extend to two wavelenghts

  z=f1(:,1); 
  ind = find(z == min(z));
  ind = ind(1);
  x1=t1(ind);
  z1=z(ind);
  p1=f1(ind,2);

  z=f2(:,1);
  ind = find(z == max(z));
  ind=ind(1);
  x2=t2(ind);
  z2=z(ind);
  p2=f2(ind,2);

  dp = p2-p1;
  l = -x1 + x2;
  L  = 4*l;
  dx = L/N;

  figure(1); plot(t1, f1(:,1), t2, f2(:,1));
      axis([x1,x2, z1, z2]); set(gca, "fontsize", 20);

%return

t =   [fliplr(t1), t2(2:end)];
z =   [flipud(f1(:,1)); f2(2:end,1)];
psi = [flipud(f1(:,2)); f2(2:end,2)];

  x = linspace(x1, x2, N/4 + 1);


figure(1); plot(1:length(t), t); set(gca, "fontsize", 20);
 

  z1   = interp1(t, z, x, 'cubic');
  s1   = interp1(t, psi, x, 'cubic');


  plot (t,z, x, z1, 'o');

%return

  x4 = [0:N-1]*dx;
  z4 = [z1, fliplr(z1(2:end-1)), z1, fliplr(z1(2:end-1) )];  
  s4 = [s1, -fliplr(s1(2:end-1))+2*p2, s1+2*dp, -fliplr(s1(2:end-1))+2*p2+dp*2];  


  delta_psi = 4*dp

  u4 = a - b./z4;

  figure(1); plot(x4,z4);
     set(gca, "fontsize", 20);
  figure(2); plot(x4,u4);
     set(gca, "fontsize", 20)

     
%return

%-- save profiles as text -- 
     
  fid=fopen([fbase, '.txt'], 'wt');
  fprintf(fid, '# Profile of travelling capillary wave from ODE\n');
  fprintf(fid, '# Parameters: a = %f (speed), b = %f,  c = %f\n', a,b,c);
  fprintf(fid, '# L = %20.10e (domain length), delta psi = %d\n\n', L, delta_psi);

  fprintf(fid, '# 1.x  2.z  3.u\n\n');
  for i=1:N
	  fprintf(fid, '%16.6e  %16.6e  %16.6e\n', x4(i), z4(i), u4(i));
  end
  fclose(fid);

return

%-- save initial conditions for 2D run --


  psi4 = (z4+1) + iunit*s4;

  psi = repmat(psi4, N, 1);

  p  = zeros(2,N,N);
  p(1,:,:) = real(psi);
  p(2,:,:) = imag(psi);

  p = reshape(p, N*N*2, 1);

  fid = fopen([fbase, '.psi.0000'], 'wb');
  fwrite(fid, p, 'double');
  fclose(fid);



end

%---------------------

function f = fdot(x,t)

  global a;
  global b;
  global c;

  f = [0; 0];

  z=x(1);


  f(1) = real( sqrt(   2*c*a*a - a*z - b*b/z  ) );
  f(2) = a * (1 - b/z);
 
end

%---------------------


