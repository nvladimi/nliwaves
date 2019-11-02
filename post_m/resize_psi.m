function resize_psi(fname_in, fname_out, N1, N2)
%
  %   resize_psi(fname_in, fname_out, N1, N2)
%
%   Input:
%
%       fname_in:   binary input file (size N1)
%       fname_out:  binary output file (size N2)
%

  fid = fopen(fname_in, 'rb');
  p = fread(fid, N1*N1*2, 'double');
  fclose(fid);

  p = reshape(p,2,N1,N1);
  p = p(1,:,:)+ i*p(2,:,:);
  p = reshape(p,N1,N1);

  fp1 = fft2(p);

  fp1(1,1)

  fp1 = fftshift(fp1);

  p1max = max(max(abs(p)));

  clear p;

  if (N2 < N1)
    ind = (N1/2 - N2/2 + 1) : (N1/2 + N2/2);
    fp2 = fp1(ind,ind);
  else 
    ind = (N2/2 - N1/2 + 1) : (N2/2 + N1/2);
    fp2 = zeros(N2,N2);
    fp2(ind,ind) = fp1;
  end

  clear fp1;

  fp2 = ifftshift(fp2);

  fp2(1,1)

  p2 = (N2/N1)**2 * ifft2(fp2);

  p2max = max(max(abs(p2)));

  p  = zeros(2,N2,N2);
  p(1,:,:) = real(p2);
  p(2,:,:) = imag(p2);

  p=reshape(p,N2*N2*2,1);

  fid = fopen(fname_out, 'wb');
  fwrite(fid, p, 'double');
  fclose(fid);

  p1max, p2max

end

