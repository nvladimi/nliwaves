
function reduce_psi
%
%   reduce_psi - see and edit the script
%

  fbase1 = 'eps1e-2_4096';
  fbase2 = 'TMP/zoom_256';

  for n=101:172

    fname_in  = [fbase1, '.psi.',  num2str(n,'%04d')];
    fname_out = [fbase2, '.psi.',  num2str(n,'%04d')];

    reduce_psi_one(fname_in, fname_out, 4096, 256);

  end

end

function reduce_psi_one(fname_in, fname_out, N, Nnew)
  % fname_in:   binary input file (size N)
  % fname_out:  binary output file (size Nnew)
  % N:          size of original data
  % Nnew:       size of reduced data

  scale=round(N/Nnew);
  ind=scale:scale:N;


  f = fopen(fname_in, 'rb');
  p = fread(f, N*N*2, 'double');
  fclose(f);

  p = reshape(p,2,N,N);

  p1 = p(:,ind,ind);
  
  N1 = length(p1);

  p1=reshape(p1,N1*N1*2,1);

  f = fopen(fname_out, 'wb');
  fwrite(f, p1, 'double');
  fclose(f);

end

