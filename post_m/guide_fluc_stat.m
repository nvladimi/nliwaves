function guide_fluc_stat
%
%  function guide_fluc_stat
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
  fbase='a400e-3';  N = 1024;
 %fbase='n8000';    N = 1024;

 %fbase='t75a12';   N = 1024;
 %fbase='t65a16';   N = 1024;
 %fbase='t40a10';   N = 1024;
 %fbase='t25a14';   N = 1024;
 %fbase='t15a20';   N = 1024;

  %--------------

  NN=N*N;


  fid = fopen([fbase,'.stat'], 'wt');

  fprintf(fid, "% 1.file  2.N  3.N0  4.|dphi|^2/N0  5.flatness\n\n");


for k=[2,50,80, 128, 130:150]
 
    num = num2str(k,'%04d');

    f = read_psi([fbase, '.psi.', num], N);

    favg = sum(sum(f))/NN;

    Nwaves0   =  abs(favg).^2;
    Nwaves    =  sum(sum(f.*conj(f)))/NN;

    s2 = sum(sum( (abs(f - favg)).^2  ))/NN;
    s4 = sum(sum( (abs(f - favg)).^4  ))/NN;
    
    flatness = s4/s2.^2;

    fprintf(fid,"%4d  %8.2f  %8.2f  %8.6f  %8.6f\n", ...
             k, Nwaves, Nwaves0, s2/Nwaves0, flatness);
    
  end

  fclose(fid);

end


%---------------------------------------------------------
