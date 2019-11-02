function psi_corrR_phase(fname_in,fname_out, N, kmax, corrfun_type)
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


  f = read_psi(fname_in, N);

  f = angle(f);

  %-- angle averaging -------------------------------------------

  kkmax = kmax*kmax;

  if (corrfun_type == 'avg')

    S=zeros(kkmax+1,1);
    C=zeros(kkmax+1,1);

    for i=0:kmax
      for j=-kmax:kmax

        kk = i*i + j*j;

        if (kk > kkmax) 
	  continue;
        end
      
        df = f - circshift(f,[-j,-i]);

        ind1=find(df > pi);
        ind2=find(df <= -pi);

        df(ind1) = df(ind1) - 2*pi;
        df(ind2) = df(ind2) + 2*pi;

        q   = df.*df;
        s   = sum(sum(q))/N/N;

        S(kk+1) = S(kk+1) + s; 
        C(kk+1) = C(kk+1) + 1;

      end
    end

    %-- writing output --

    fid = fopen(fname_out,'wt');

    fprintf(fid,'%% Angle-averaged correlation functions of phase \n');
    fprintf(fid,'%% computed by \"psi_corrR_phase.m\" from \"%s\".\n', fname_in);
    fprintf(fid,'%%\n% 1.index  2.corrfun\n\n');

    for kk=0:kkmax
      if (C(kk+1) ~= 0)
        fprintf(fid,'%8.4f  %16.8e\n', sqrt(kk), S(kk+1)/C(kk+1));
      end
    end

    fclose(fid);


  %-- anisotropy ------------------------------------------------

  else if (corrfun_type == 'dir')

    S=zeros(kmax+1,4);

    s1 =  [1,0];  % x-direction
    s2 =  [0,1];  % y-direction
    s3 =  [1,1];  % xy-direction
    s4 = [-1,1];  % yx-direction

    for k=0:kmax

      %-- x-direction 

      i=s1(1)*k;
      j=s1(2)*k;

      df = f - circshift(f,[-j,-i]);

      ind1=find(df > pi);
      ind2=find(df <= -pi);

      df(ind1) = df(ind1) - 2*pi;
      df(ind2) = df(ind2) + 2*pi;

      q   = df.*df;
      s   = sum(sum(q))/N/N;

      S(k+1,1) = s; 
 

      %-- y-direction 

      i=s2(1)*k;
      j=s2(2)*k;

      df = f - circshift(f,[-j,-i]);

      ind1=find(df > pi);
      ind2=find(df <= -pi);

      df(ind1) = df(ind1) - 2*pi;
      df(ind2) = df(ind2) + 2*pi;

      q   = df.*df;
      s   = sum(sum(q))/N/N;

      S(k+1,2) = s; 


      %-- xy-direction 

      i=s3(1)*k;
      j=s3(2)*k;

      df = f - circshift(f,[-j,-i]);

      ind1=find(df > pi);
      ind2=find(df <= -pi);

      df(ind1) = df(ind1) - 2*pi;
      df(ind2) = df(ind2) + 2*pi;

      q   = df.*df;
      s   = sum(sum(q))/N/N;

      S(k+1,3) = s; 


      %-- yx-direction 

      i=s4(1)*k;
      j=s4(2)*k;

      df = f - circshift(f,[-j,-i]);

      ind1=find(df > pi);
      ind2=find(df <= -pi);

      df(ind1) = df(ind1) - 2*pi;
      df(ind2) = df(ind2) + 2*pi;

      q   = df.*df;
      s   = sum(sum(q))/N/N;

      S(k+1,4) = s; 

    end


    %-- writing output --

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

  else 
     error('unknown corrfun_type')
  end

















end

