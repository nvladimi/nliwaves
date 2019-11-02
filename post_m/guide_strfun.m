function s2=guide_strfun(f,N)
%
%  function guide_strfun
%
%  Paramemters are edited inside function.
%
%  Reads *.psi files and creates images of magnitude (deviation 
%  from average), phase, and spectrum n_k without dealiased modes.
%  Fixed colomaps.  Computation of qflux is commented out.

%----------------------------------------------------------

% for debuggung and testing: [s2,s4]=guide_strfun(f,N)
  s2=strfun2(f,N); return; s4=strfun4(f,N); return;


 %fbase='a4';       N =  256;
 %fbase='a064b';    N =  256;
 %fbase='a800e-3';  N =  512;
 %fbase='a160e-2';  N =  512;
 %fbase='a400e-3';  N = 1024;

 %fbase='t75a12';   N = 1024;
 %fbase='t65a16';   N = 1024;
 %fbase='t40a10';   N = 1024;
 %fbase='t25a14';   N = 1024;
  fbase='t15a20';   N = 1024;
 %fbase='t01a25';   N = 1024;


  %--------------


  s2=zeros(N,N);
  s4=zeros(N,N);

  n = 0;

  for k=51 %:150
 
    num = num2str(k,'%04d');

    %[fbase, '.psi.', num]

    f = read_psi([fbase, '.psi.', num], N);

    s2 = s2 + strfun2(f,N);
    s4 = s4 + strfun4(f,N);

    n  = n + 1;

  end


  s2 = s2/n;
  s4 = s4/n;

  fid=fopen([fbase, '.sf2.avg'], 'wb');
  fwrite(fid, s2, 'double');
  fclose(fid);

  fid=fopen([fbase, '.sf4.avg'], 'wb');
  fwrite(fid, s4, 'double');
  fclose(fid);

end



%---------------------------------------------------------

function s = strfun2(psi, N)

    c1 = sum( sum(psi.*conj(psi)) );

    f = psi;
    g = f;
    c2 = ifftn( conj(fftn(f)) .* fftn(g) );

    s = real( 2*c1 - 2*c2 );

    s = s/N.^2;


end

%---------------------------------------------------------

function s = strfun4(psi, N)


    c1 = sum(sum( (psi.*conj(psi)).^2  ));

    f = psi.*psi;
    g = f;
    c2 = ifftn( conj(fftn(f)) .* fftn(g) );

    f = psi.*conj(psi);
    g = f;
    c3 = ifftn( conj(fftn(f)) .* fftn(g) );

    f = psi.*psi.*conj(psi);
    g = psi;
    c4 = ifftn( conj(fftn(f)) .* fftn(g) );

    f = psi;
    g = psi.*psi.*conj(psi);
    c5 = ifftn( conj(fftn(f)) .* fftn(g) );

    s = real( 2*c1 + 2*c2  + 4*c3 - 4*c4 - 4*c5 );

    s = s/N.^2;


end

%---------------------------------------------------------
