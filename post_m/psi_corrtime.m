
function tcorr = psi_corrtime(fname_in, fname_out, N, M, dt)
%
%   psi_corrtime(fname_in, fname_out, N, M, dt)
%
%   Compute time correlation function of Psi.
%
%   Input:
%        fname_in       binary input file with psi, NxNxM
%        fname_out      text output file
%        N              number of spacial points is NxN
%        M              number of time layers  
%        dt             time interval between layers


  %-- read data in and shape as NxNxM array --  

  fid = fopen(fname_in, 'rb');
  p   = fread(fid, N*N*M*2, 'double');
  fclose(fid);

  p = reshape(p,2,N,N,M);

  u = p(1,:,:,:);
  v = p(2,:,:,:);

  u = reshape(u,N,N,M);
  v = reshape(v,N,N,M);

  p = u+i*v;

  %-- compute correlation function for each point --

  fp = fft(p,M,3);

  fp = fp.*conj(fp);

  p  = ifft(fp,M,3)/M;
  p0 = real(p(:,:,1));

  for k=1:M
     p(:,:,k)=p(:,:,k)./p0;
  end


  %-- average over points --

  s = sum(sum(p))/N/N;
  s = reshape(s,[1,M]);
  

  %-- correlation time -- 

  i1 = max(find(real(s(1:M/2))>0.5));
  s1 = real(s(i1));
  s2 = real(s(i1+1));
    
  tcorr = (i1 + (0.5-s1)/(s2-s1))*dt;


  %-- write to file --

  fid = fopen(fname_out,'wt');

  fprintf(fid,'%% Time correlation function of Psi \n');
  fprintf(fid,'%% computed by \"psi_corrtime.m\" from \"%s\".\n', fname_in);
  fprintf(fid,'%% correlation time T = %12.4e\n', tcorr);
  fprintf(fid,'%%\n%% 1.index  2.time  3.Re(corrfun)  4.Im(corrfun)\n\n');

  t=M*dt;

  for k=1:M/2 
    fprintf(fid,'%6d  %10.6f  %16.8e  %16.8e\n', ...
            k, (k-1)*dt, real(s(k)),  imag(s(k)) );
  end

  fprintf(fid,'\n\n');

  for k=M/2+1:M; 
    fprintf(fid,'%6d  %10.6f  %16.8e  %16.8e\n', ...
            k, (k-1)*dt-t, real(s(k)),  imag(s(k)) );
  end


  fclose(fid);

end

