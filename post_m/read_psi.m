function f = read_psi(fname, N)
%
%   f = read_psi(fname, N)
%
%   reads file, returns complex array NxN
%
%   Input:
%        fname   file name to read
%        N       domain size

  %fname = '../code/test/gauss.psi.0000';
  %N = 256;

  f = fopen(fname, 'rb');
  p = fread(f, N*N*2, 'double');
  fclose(f);

  p = reshape(p,2,N,N);

  u = p(1,:,:);
  v = p(2,:,:);

  u = (reshape(u,N,N))';
  v = (reshape(v,N,N))';

  f = u + i*v;

  %imagesc(u);

end

