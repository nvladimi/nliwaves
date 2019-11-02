function guide_slice_angular
%
%  function guide_slice_angular
%
%  Paramemters are edited inside function.
%
%  Read *nk file and writes text file with angular dependence of nk
%  along several circles.

%----------------------------------------------------------

  global F;
  global N;

  %infile='a064b.psi.0090';  N = 256;

  %infile='s115.nk.avg';  N = 512;    % t = 3000
  %infile='s17.nk.avg';   N = 512;    % t = 6000
  %infile='s16.nk.avg';   N = 512;    % t = 7500

  %infile='s20.nk.avg';   N = 1024;   % t = 1500
  %infile='s27.nk.avg';   N = 1024;   % t = 2100

  %infile='t75a12.nk.avg';   N = 1024;   % t = 7500
   infile='t65a16.nk.avg';   N = 1024;   % t = 6500
  %infile='t40a10.nk.avg';   N = 1024;   % t = 4000
  %infile='t25a14.nk.avg';   N = 1024;   % t = 2500
  %infile='t15a20.nk.avg';   N = 1024;   % t = 1500
  %infile='t01a25.nk.avg';   N = 1024;   % t =  100

  amax=180;  % number of points per circle
  rmax=8;    % number of circular slices

  outfile = [infile, '.aslices']


%-----------------------

  F = read_fft(infile);
  %F = make_fft(infile);


%-----------------------

  a = 2*pi*(0:amax-1)/amax; 
  r = (1:rmax)/rmax*N/8;

  f = nk_slice(r, a);

  fmax = max(f);
  fsum = sum(f);

  fid=fopen(outfile, 'wt');

  fprintf(fid, '#1.a_2.r=  ');
  fprintf(fid, '%f  ', r);
  fprintf(fid, '\n');

  fprintf(fid, '#max =     ');
  fprintf(fid, '%10.4e  ', fmax);
  fprintf(fid, '\n\n');

  fprintf(fid, '#sum =     ');
  fprintf(fid, '%10.4e  ', fsum);
  fprintf(fid, '\n\n');

  for i=1:amax
    fprintf(fid, '  %10.4e', a(i) );
    fprintf(fid, '  %10.4e', f(i,:) );
    fprintf(fid, '\n');
  end

  fclose(fid);

end

%---------------------------------------------------------


function zi = nk_slice(r, a)

  global F;
  global N;

  [rr,aa] = meshgrid(r,a);

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
  s = fftshift(reshape(f,N,N));

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

