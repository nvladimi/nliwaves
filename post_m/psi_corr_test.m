function [Npart, r] = psi_corr_test %(fname_in,fname_out, N, kmax, corrfun_type)
%
%   psi_corrR_phase(fname_in,fname_out, N, kmax, corrfun_type)
%
%   Compute correlation function of phase. Corelation function 
%   is either angle ageraged, or computed in x,y, xy, or yx directions.
%
%   Input:
%        fname_in       binary input file with psi
%        fname_out      text output file
%        N              size of simulation
%        kmax           extend of the correlation function in points
%        corrfun_type   'avg' or 'dir'

  fname_in = 'f1_a0400.psi.0040';
  N = 1024;
  L = 2*pi;
  
  kmax = 6;
  

  %------------------------------------

  f  =  read_psi(fname_in, N);
  F0 =  zeros(N,N);
  F2 =  zeros(N,N);
  dx = L/N;

  %imagesc(abs(F)); 

  f0 = conj(f);

  kkmax = kmax*kmax;

  for i=-kmax:kmax
     for j=-kmax:kmax

        kk = i*i + j*j;

        if (kk > kkmax) 
	  continue;
        end
      
        dF = f0.*circshift(f,[-j,-i]); 
 
        F0 = F0 + dF; 
        F2 = F2 + dF*kk; 

     end
  end

  F0 = F0*dx^2;
  F2 = F2*dx^4;

  f = 0.5*abs(F0)/11.7;
  r = sqrt(abs(F2)./abs(F0));
 
  f1=circshift(f,[0,1]); f2=circshift(f,[0,-1]); f3=circshift(f,[1,0]); f4=circshift(f,[-1,0]);
  ind=find(f>f1 & f>f2 & f>f3 & f>f4 & f > 1); 

  Npart = f(ind); 
  r     = r(ind);

 


    %-- writing output --

    return;

    fid = fopen(fname_out,'wt');
 
    fprintf(fid,'%% Directional correlation functions of phase\n');
    fprintf(fid,'%% computed by \"psi_corrR_phase.m\" from \"%s\".\n', fname_in);
    fprintf(fid,'%%\n%% 1.index  2.x-dir  3.y-dir  4.xy-dir  5.yx-dir\n\n');

    for k=0:kmax
      fprintf(fid,'%6d', k);
      fprintf(fid,' %16.8e', S(k+1, :));
      fprintf(fid,'\n');
    end

    fclose(fid);

  %-------------------------------------------------------------------



end

