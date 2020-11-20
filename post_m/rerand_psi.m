function rerand_psi(fname_in, fname_out, N, seed, randomness)
%
%   rerand_psi(fname_in, fname_out, N, seed, randomness)
%
%   Input:
%
%       fname_in:   binary input file (size N)
%       fname_out:  binary output file (size N)
%       N:          number of gridpoints in each direction
%       seed:       master seed (integer)
%       randomness: [0,1)
%

  fid = fopen(fname_in, 'rb');
  p = fread(fid, N*N*2, 'double');
  fclose(fid);

  p = reshape(p,2,N,N);
  p = p(1,:,:)+ 1i*p(2,:,:);
  p = reshape(p,N,N);


  rand ("twister", seed);

  r = (rand(N,N)*2 - 1)*pi*randomness;

  fp = fft2(p);
  fp = fp .* exp(1i*r);
  pnew = ifft2(fp);

  p  = zeros(2,N,N);
  p(1,:,:) = real(pnew);
  p(2,:,:) = imag(pnew);

  p=reshape(p,N*N*2,1);

  fid = fopen(fname_out, 'wb');
  fwrite(fid, p, 'double');
  fclose(fid);


end

