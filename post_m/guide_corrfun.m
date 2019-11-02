function guide_corrfun
%
%  function guide_images_quick
%
%  Paramemters are edited inside function.
%
%  Reads *.psi files and creates images of magnitude (deviation 
%  from average), phase, and spectrum n_k without dealiased modes.
%  Fixed colomaps.  Computation of qflux is commented out.

%----------------------------------------------------------

 %fbase='a4';       N =  256;
 %fbase='a064b';     N =  256;
 %fbase='a800e-3';  N =  512;
 %fbase='a160e-2';  N =  512;
  fbase='a400e-3';  N = 1024;  amax=8; a0= 0.25;

  %infile='a064b';  N = 256;  amax=8; a0=-0.3+pi/4; 

  %--------------

  rmax=N/2;  % number of points per slices

  a = a0 + 2*pi*(0:amax-1)/amax; 
  r = (0:rmax-1)/rmax*N/2;

  %--------------

  for k=150
 
    num = num2str(k,'%04d');

    f = read_psi([fbase, '.psi.', num], N);

    f = abs(f) - sum(sum(abs(f)))/(N*N);

    %-- <psi*conj(psi)> --

    s = fftn(f) .* conj(fftn(f));
    s = ifftn(s)/(N*N);

    s = real(s);
    s = s/s(1,1);

    s = fftshift(s);

    %one_image(s, -0.5, 1, [fbase, '.cf.', num, '.png']);


    %-----------------------

    outfile = [fbase, '.cf.', num, '.rslices']

    f = nk_rslice(s, N, a, r);

    fid=fopen(outfile, 'wt');

    fprintf(fid, '#1.r_2.a=  ');
    fprintf(fid, '%f  ', a);
    fprintf(fid, '\n');

    for i=1:rmax
      fprintf(fid, '  %10.4e', r(i) );
      fprintf(fid, '  %10.4e', f(i,:) );
      fprintf(fid, '\n');
    end

    fclose(fid);

  end

end


%---------------------------------------------------------

function s = read_fft(fname, N)

   fid=fopen(fname, "rb");  
   f = fread(fid, N*N, 'double');  
   fclose(fid);  
   s = fftshift(reshape(f,N,N));
   %s = circshift(reshape(f,N,N), [N/4,N/4]);

end


%---------------------------------------------------------

%---------------------------------------------------------

function zi = nk_rslice(F, N, a, r)

  [aa,rr] = meshgrid(a,r);

  xi = rr.*cos(aa) + N/2 + 1;
  yi = rr.*sin(aa) + N/2 + 1;

  zi = interp2(F,xi,yi, 'linear');

end 

%---------------------------------------------------------
