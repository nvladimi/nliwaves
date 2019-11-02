function guide_nk_avg
%
%  function guide_nk_avg
%
%  Paramemters are edited inside function.
%
%  Reads *.psi files and creates images of magnitude (deviation 
%  from average), phase, and spectrum n_k without dealiased modes.
%  Fixed colomaps.  Computation of qflux is commented out.

%----------------------------------------------------------

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
  fbase='t25a14';   N = 1024;
 %fbase='t15a20';   N = 1024;

 %fbase='t01a25';   N = 1024;
  %--------------

  NN = N*N;
  nk = zeros(N,N);
  nf = 0;


  for k=1:100
 
    num = num2str(k,'%04d');

    f = read_psi([fbase, '.psi.', num], N);

    q=fft2(f)/NN; q=q.*conj(q);  

    nk = nk + real(q);

    nf = nf+1;
    
  end

  nk = nk/nf;

  fname=[fbase, '.nk.avgm'];
  fid=fopen(fname, 'wb');
  fwrite(fid, nk, 'double');
  fclose(fid);  
 
end


%---------------------------------------------------------
