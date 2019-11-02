
function guide_nkflow
%
%  function guide_images_quick
%
%  Parameters are edited inside.
%
%  Reads *.nk files, computes forcing, source term, and potential
%  for visualizing q-flux.  Data is saved in vtk format for paraview.
%  
%  Subroutine nk_shells prints on the screen source term integrated
%  over pumping and damping shells.  
%  


 %fbase='a064';      N = 128;

 %fbase='a160e-2s';  N = 256;   % t=1000
 %fbase='s115';      N = 256;   % t=3000
 %fbase='s17';       N = 256;   % t=6000
 %fbase='s20';       N = 512;   % t=1500
 %fbase='s27';       N = 512;   % t=2100

 %fbase='t75a12';   N = 1024/2;
 %fbase='t65a16';   N = 1024/2;
 %fbase='t40a10';   N = 1024/2;
 %fbase='t25a14';   N = 1024/2;
  fbase='t15a20';   N = 1024/2;


  %-- create charges ---

  k = 0;
  %num = num2str(k,'%04d');
  num='avgm';

  fname = [fbase, '.nk.', num];

  nk = read_fft(fname, N);
  fk = guide_forcing(2, N);

  qk = nk.*fk;

  %--debugging--
  %
  %fmin=min(min(q)), fmax=max(max(q)), 
  %one_image(f,  -5, 5, 'force.png');
  %one_image(q,  -0.1, 0.1, 'source.png');

  nk_shells(qk,N);  

  qtot = sum(sum(qk));

  qk(N/2+1,N/2+1) = -qtot;
  %qk(N/2+1,N/2+1) = 0;    % debugging?

  %return;

%-- compute potential (expensive) -- 

  pk = zeros(N,N);

  [x, y] = meshgrid(1:N, 1:N);
  for i= 1:N
    for j=1:N
       s = qk ./ sqrt((x-i).^2 + (y-j).^2);
       s(j,i) = 0;
       pk(j,i) = sum (sum(s)); 
    end
  end


  % solve poisson equation spectrally
  %
  %rho = fft2(qk);
  %kx = [0:N/2,-N/2+1:-1];
  %[kkx, kky] = meshgrid(kx, kx);
  %kk = kkx.*kkx + kky.*kky;
  %pk = rho./kk;
  %pk(1,1) = 0;
  %pk = real(ifft2(phi));


  fmin=min(min(pk)); fmax=max(max(pk));
  one_image(pk,  fmin, fmax, 'potential_1.png');



  %-- save in files --

  vtk_any(nk, 'nk', [fbase, '_nk.vtk']);
  vtk_any(qk, 'qk', [fbase, '_qk.vtk']);
  vtk_any(pk, 'pk', [fbase, '_pk.vtk']);


end

%--------------------------------------------------

function s = read_fft(fname, N)

   ind = N/2+1 : 3*N/2;

   fid=fopen(fname, "rb");  
   f = fread(fid, 4*N*N, 'double');  
   fclose(fid);  
   s = fftshift(reshape(f, 2*N, 2*N));
   s = s(ind,ind);

end




%--------------------------------------------------

  function nk_shells(f,N)
% N is number of dealiased modes.

  s = N/128;

  n  = 128*s;
  kl =  28*s;
  kr =  32*s;

 
  %-- grid --

  k = [-n/2:n/2-1]; 

  [kx, ky] = meshgrid (k, k);

  k = sqrt(kx.^2 + ky.^2);

  %-- find total --

  q3 = sum(sum(f));

  %-- erase damping shell --

  ind = find (k>=kr-1);
  f(ind) = 0;

  q2 = sum(sum(f));


  %-- erase pumping shell --

  ind = find (k>kl);
  f(ind) = 0;

  q1 = sum(sum(f));

  %---------------------

  disp("   total  inner  pumping  dumping");

  disp([q3, q1, q2-q1, q3-q2]); 

end

%--------------------------------------------------



