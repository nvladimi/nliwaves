function guide_images_all
%
%  function guide_images_all
%
%  Paramemters are edited inside function.
%
%  Reads *psi, *nk, and *sfN files for condensate and fluctuations.
%  From *psi files, genarates images of magnitude in  adjustable colomap,
%  and phase, in fixed colomap.
%
%  From *nk files generates images of log(nk) in fixed colomap. 
%
%  Images of correlation functions and flatness are commented out.

%----------------------------------------------------------

 %fbase='a064';     N =  256;
 %fbase='tmp';      N =  256;
 %fbase='a160e-2s'; N =  512;
 %fbase='t6000';    N =  512;
 %fbase='s115';     N =  512;   % t=3000
 %fbase='s17';      N =  512;   % t=6000
 %fbase='s16';      N =  512;   % t=7500
 %fbase='a400e-3';  N = 1024;
 %fbase='s2100';    N = 1024;
 %fbase='s20';      N = 1024;   % t=1500
 %fbase='s27';      N = 1024;   % t=2100

 %fbase='t75a12';   N = 1024;
 %fbase='t65a16';   N = 1024;
 %fbase='t40a10';   N = 1024;
 %fbase='t25a14';   N = 1024;
  fbase='t15a20';   N = 1024;



  NN=N*N;


  nkmin  = log(1.e-16);        % spectrum range for log(n_k)
  nkmax  = log(1.e2);          %


for k=0 %140:10:100
 
    %num = num2str(k,'%04d');
    num='avgm';



    %-- n_k spectra --

    f = read_fft([fbase, '.nk.', num], N);

    one_image(log(f), nkmin, nkmax, [fbase, '.nk.', num, '.png']);

    continue

    f = read_fft([fbase, '_cnd.nk.', num], N);
    one_image(log(f), nkmin, nkmax, [fbase, '_cnd.nk.', num, '.png']);

    f = read_fft([fbase, '_flc.nk.', num], N);
    one_image(log(f), nkmin, nkmax, [fbase, '_flc.nk.', num, '.png']);

    continue


    %-- amplitudes and phases --

    %f = read_psi([fbase, '.psi.', num], N);
    %fmin=min(min(abs(f))); fmax=max(max(abs(f)));
    %one_image(abs(f),  fmin, fmax, [fbase, '.abs.', num, '.png']);
    %one_image(arg(f), -pi, pi, [fbase, '.arg.', num, '.png']);
    %disp([fmin, fmax])
   
    %f = read_psi([fbase, '_cnd.psi.', num], N);
    %fmin=min(min(abs(f))); fmax=max(max(abs(f)));
    %one_image(abs(f),  fmin, fmax, [fbase, '_cnd.abs.', num, '.png']);
    %one_image(arg(f), -pi,pi, [fbase, '_cnd.arg.', num, '.png']);
    %disp([fmin, fmax])

    %f = read_psi([fbase, '_flc.psi.', num], N);
    %fmax=max(max(abs(f)));
    %one_image(abs(f),  0,  fmax, [fbase, '_flc.abs.', num, '.png']);
    %one_image(arg(f), -pi, pi, [fbase, '_flc.arg.', num, '.png']);
    %disp(fmax)

    %continue


    %-- correlation functions and flatness --


    s2 = read_fft([fbase, '_cnd.sf2.', num], N);
    one_image(s2, 0, 0.08,  [fbase, '_cnd.sf2.', num, '.png']);
    %disp([min(min(s2)), max(max(s2))])

    s4 = read_fft([fbase, '_cnd.sf4.', num], N);
    one_image(s4, 0, 0.015, [fbase, '_cnd.sf4.', num, '.png']);
    %disp([min(min(s4)), max(max(s4))])

    f  = s4./s2.^2;
    one_image(f,  2, 3,     [fbase, '_cnd.flt.', num, '.png']);
    %disp([min(min(f)), max(max(f))])

    s2 = read_fft([fbase, '_flc.sf2.', num], N);
    one_image(s2, 0, 60,    [fbase, '_flc.sf2.', num, '.png']);
    %disp([min(min(s2)), max(max(s2))])

    s4 = read_fft([fbase, '_flc.sf4.', num], N);
    one_image(s4, 0, 10000,  [fbase, '_flc.sf4.', num, '.png']);
    %disp([min(min(s4)), max(max(s4))])

    f  = s4./s2.^2;
    one_image(f,  2, 3,     [fbase, '_flc.flt.', num, '.png']);
    %disp([min(min(f)), max(max(f))])
 
  end

end


%---------------------------------------------------------

function s = read_fft(fname, N)
   ink = N/4+1:3*N/4;

   fid=fopen(fname, "rb");  
   f = fread(fid, N*N, 'double');  
   fclose(fid);  
   s = fftshift(reshape(f,N,N));
   s = s(ink,ink);

   %s = circshift(reshape(f,N,N), [N/4,N/4]);

end


%---------------------------------------------------------

