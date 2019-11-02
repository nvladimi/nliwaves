function psi_corrR_conv(fname_in,fname_out, N, kmax, corrfun_type)
%
%   psi_corrR_conv(fname_in,fname_out, N, kmax, corrfun_type)
%
%   Compute correlation function of two quantities using convoluton. 
%   Corelation function is either angle ageraged, or computed in 
%   x,y, xy, or yx directions.
%
%   Input:
%        fname_in       binary input file with psi
%        fname_out      text output file
%        N              size of simulation
%        kmax           extend of the correlation function in points
%        corrfun_type   'avg' or 'dir'


  f = read_psi(fname_in, N);

%-- build two quantities to correlate: f and g  ---------------

  %g = conj(f);

  conv_name = '<psi*conj(psi)>';


%-- compute convolution ---------------------------------------

  f = fftn(f) .* conj(fftn(f));

  f = ifftn(f)/(N*N);

  f = fftshift(f);
  
  i0=N/2+1;
  j0=N/2+1;

  if (kmax >= N/2)
     kmax = N/2-1;
  end

  kkmax = kmax*kmax;

%-- angle averaging -------------------------------------------

  if (corrfun_type == 'avg')

    S=zeros(kkmax+1,1);
    C=zeros(kkmax+1,1);

    for j=-kmax:kmax
      for i=-kmax:kmax

        kk = i*i + j*j;

        if (kk <= kkmax) 
          S(kk+1) = S(kk+1) + f(i0+i,j0+j); 
          C(kk+1) = C(kk+1) + 1;
        end

      end
    end

    %-- writing output --

    fid = fopen(fname_out,'wt');

    fprintf(fid,'%% Angle-averaged correlation functions, %s\n', conv_name);
    fprintf(fid,'%% computed by \"psi_corrR_conv.m\" from \"%s\".\n', fname_in);
    fprintf(fid,'%%\n%% 1.index  2.Re(corrfun)  3.Im(corrfun)\n\n');

    for kk=0:kkmax
      if (C(kk+1) ~= 0)
        fprintf(fid,'%8.4f  %16.8e  %16.8e\n', ...
		sqrt(kk), real(S(kk+1))/C(kk+1), imag(S(kk+1))/C(kk+1) );
      end
    end

    fclose(fid);


%-- anisotropy ------------------------------------------------

  else if (corrfun_type == 'dir')

    S=zeros(kmax+1,4);

    S(1,:) = f(i0,j0);

    for k=1:kmax

      %-- x-direction 

      S(k+1,1) = ( f(i0-k,j0) + f(i0+k, j0) )/2 ; 

      %-- y-direction 

      S(k+1,2) = ( f(i0,j0-k) + f(i0, j0+k) )/2 ; 

      %-- x-direction 

      S(k+1,3) = ( f(i0-k,j0-k) + f(i0+k, j0+k) )/2 ; 

      %-- x-direction 

      S(k+1,4) = ( f(i0-k,j0+k) + f(i0+k, j0-k) )/2 ; 

    end


    %-- writing output --

    fid = fopen(fname_out,'wt');
 
    fprintf(fid,'%% Directional correlation functions, %s\n', conv_name);
    fprintf(fid,'%% computed by \"psi_corrR_gradphase.m\" from \"%s\".\n', fname_in);
    fprintf(fid,'%%\n%% 1.index\n');
    fprintf(fid,'%% 2.x-dir(Re)  3.y-dir(Re)  4.xy-dir(Re)  5.yx-dir(Re)\n');
    fprintf(fid,'%% 6.x-dir(Re)  7.y-dir(Re)  8.xy-dir(Re)  9.yx-dir(Re)\n\n');
 
    for k=0:kmax
      fprintf(fid,'%6d', k);
      fprintf(fid,' %16.8e', real(S(k+1, :)) );
      fprintf(fid,' %16.8e', imag(S(k+1, :)) );
      fprintf(fid,'\n');
    end

    fclose(fid);

  %-------------------------------------------------------------------

  else 
     error('unknown corrfun_type')
  end

end


