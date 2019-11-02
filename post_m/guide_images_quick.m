function guide_images_quick
%
%  function guide_images_quick
%
%  Paramemters are edited inside function.
%
%  Reads *.psi files and creates images of magnitude (deviation 
%  from average), phase, and spectrum n_k without dealiased modes.
%  Fixed colomaps.  Computation of qflux is commented out.

%----------------------------------------------------------

%fbase='re_f3/n3600';  N=256;
fbase='f0001';  N=256;

%fbase='n0800';  N=256;

 %fbase='a1';       N =  256;
 %fbase='a064c';    N =  256;
 %fbase='n6000';    N =  256;
 %fbase='a800e-3';  N =  512;
 %fbase='a160e-2';  N =  512;
 %fbase='a400e-3';  N = 1024;
 %fbase='n8000';    N = 1024;

 %fbase='t75a12';   N = 1024;
 %fbase='t65a16';   N = 1024;
 %fbase='t40a10';   N = 1024;
 %fbase='t25a14';   N = 1024;
 %fbase='t15a20';   N = 1024;

  %--------------

  nkmin  = log(1.e-16);        % spectrum range for log(n_k)
  nkmax  = log(1.e2);         %

  %nkmin  = log(1.e-12)
  %nkmax  = log(2.e2);         %

  %nkmin  = log(1.e-32);        % spectrum range for log(n_k)
  %nkmax  = log(1.e-4);         %

  df     = 3;          % |psi| range = avg +/- df


  %--------------

  NN=N*N;
  ink = N/4+1:3*N/4;

for k=0:10  %[2,30,40,42,50,80,130, 142:2:148]
 
    num = num2str(k,'%04d');

    f = read_psi([fbase, '.psi.', num], N);
    %p = phase_diff(f);

    absf = abs(f);
    avgf = sum(sum(absf))/NN;
    %one_image(absf-avgf,  -df, df, [fbase, '.abs.', num, '.png']);
    %one_image(arg(f), -pi, pi, [fbase, '.arg.', num, '.png']);
    %imagesc(absf-avgf,  [-1, 1]);  axis('equal', 'tics', 'off');


    q=fft2(f)/NN;
    q=q.*conj(q); 
    nk=q(1,1);
    q=fftshift(q);  q=q(ink,ink);
    %one_image(log(q), nkmin, nkmax, [fbase, '.nk.', num, '.png']);
    %imagesc(log(q), [nkmin, nkmax]); axis('equal', 'tics', 'off');

   k = -N/4+1:N/4;
  loglog(k, q(N/4+1,:), k, 1./k )

    %imagesc(p, [0, 2*pi]); axis('equal', 'tics', 'off');

    pause(2.0);

    printf("%4d: nk0 = %8.2f   avg = %8.2f  diff= %8.2f %8.2f\n", ...
         k, nk, avgf, min(min(absf))-avgf, max(max(absf))-avgf);
    
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

function th =  phase_diff(f)

  N=length(f);

  NN=N*N; ink = N/4+1:3*N/4;

  f=fft2(f)/NN;  a0=arg(f(1,1));  f=fftshift(f);  f=f(ink,ink); f=f(2:end, 2:end); 
   n=(abs(f)).^2; a1=arg(f); a2=flipud(fliplr(a1)); th=2*a0 - a1 - a2;

  ind=find(th<0); th(ind)=th(ind)+2*pi;
  ind=find(th<0); th(ind)=th(ind)+2*pi;

  imagesc(th, [0,2*pi])

end

%-------------------------------------------------------------
