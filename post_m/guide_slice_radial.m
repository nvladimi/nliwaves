function guide_slice_radial
%
%  function guide_slice_radial
%
%  Paramemters are edited inside function.
%
%  Read *nk file and writes text file with radial dependence of nk
%  in several directions.

%----------------------------------------------------------

  global F;
  global N;

  % amax - number of angles
  % a0 - first angle

  %infile='a160e-2s_cnd.sf4.avg';  N = 512;  amax=4; a0= 0.07;   % t = 1000
  %infile='s115_flc.sf4.avg';      N = 512;  amax=6; a0= 0.05;   % t = 3000
  %infile='s17.nk.avg';            N = 512;  amax=6; a0= 0.00;   % t = 6000
  %infile='s16_cnd.sf4.avg';       N = 512;  amax=8; a0= 0.60+pi/4; % t = 7500
  %infile='s27.nk.avg';           N = 1024;  amax=6; a0= 0.05;   % t = 2100
  %infile='a064c.psi.0080';        N = 256;  amax=8; a0=-0.24+pi/4; 

  %infile='t75a12.sf4.avg';       N = 1024;  amax=8; a0= 0.25;   % t = 7500
  %infile='t65a16.nk.avg';        N = 1024;  amax=6; a0= 1.75;   % t = 6500
   infile='t25a14_cnd.nk.avg';       N = 1024;  amax=6; a0= 1.20;   % t = 2500
  %infile='t15a20.sf4.avg';       N = 1024;  amax=4; a0= 2.67;   % t = 1500


 %fbase='t65a16';   N = 1024;
 %fbase='t40a10';   N = 1024;
 %fbase='t15a20';   N = 1024;


  rmax=512;  % number of points per slices

  outfile = [infile, '.rslices']

%-----------------------

  F = read_fft(infile);
  %F = make_fft(infile);

  a = a0 + 2*pi*(0:amax-1)/amax; 
  r = (1:rmax)/rmax*N/2;

  f = nk_rslice(a, r);

  fmax = max(f);

  fid=fopen(outfile, 'wt');

  fprintf(fid, '#1.r_2.a=  ');
  fprintf(fid, '%f  ', a);
  fprintf(fid, '\n');
  fprintf(fid, '#max =     ');
  fprintf(fid, '%10.4e  ', fmax);
  fprintf(fid, '\n\n');

  for i=1:rmax
    fprintf(fid, '  %10.4e', r(i) );
    fprintf(fid, '  %10.4e', f(i,:) );
    fprintf(fid, '\n');
  end

  fclose(fid);

end

%---------------------------------------------------------


function zi = nk_rslice(a,r)

  global F;
  global N;

  [aa,rr] = meshgrid(a,r);

  xi = rr.*cos(aa) + N/2 + 1;
  yi = rr.*sin(aa) + N/2 + 1;

  zi = interp2(F,xi,yi, 'linear');

end 

%---------------------------------------------------------

function s = read_fft(fname)

  global N;

  fid=fopen(fname, "rb");  
  f = fread(fid, N*N, 'double');  
  fclose(fid);  
  s = reshape(f,N,N);
  s = fftshift(s);

end

%---------------------------------------------------------

function q = make_fft(fname)

  global N;

  q = read_psi(fname, N);
  q = fft2(q)/N/N; 
  q = q.*conj(q);  
  q = fftshift(q);

end

%---------------------------------------------------------

