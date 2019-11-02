
function w = weights_fft(N)

  % w = zeros(N,1);
  % for j= -N:N-1
  %  for i= -N:N-1 
  %     kk = i*i + j*j;
  %     k = floor(sqrt(kk) + 0.5);
  %     if (k < N)
  %          w(k+1) = w(k+1) + 1;
  %     end
  %   end
  % end

  [kx,ky]=meshgrid(-N:N-1, -N:N-1);
  k = floor(sqrt(kx.*kx + ky.*ky) + 0.5);
  k = reshape(k,4*N*N);
  w = hist(k,0:N);
  w = w(1:N)';

end

