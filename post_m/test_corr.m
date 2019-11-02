function test_corr

  N    = 128;
  kmax = N/2;

  NN = N*N;

  r = (0:N-1)*2*pi/N;
  [x, y] = meshgrid (r, r);

  f = i*cos(x).*sin(4*y + x);

  f = circshift(f,[-30,-50]);

  figure(1); imagesc(flipud(abs(f)));

  [s2,s4]=guide_strfun(f,N);

  figure(2); imagesc(flipud(s4)); max(max(s4))

  %return

 %-----------------------


  S=zeros(2*kmax,2*kmax);
  
  o = kmax+1;

  for j=-kmax:kmax-1
    for i=-kmax:kmax-1

        df = f - circshift(f,[-j,-i]);
        s   = (df.*conj(df)).^2;

        S(j+o,i+o) = sum(sum(s))/NN;

    end
  end

  S=fftshift(S);
  S(1,1)

  figure(3); imagesc(flipud(S)); max(max(S))

end

