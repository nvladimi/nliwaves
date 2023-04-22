%------------------------------------------------------
function S = giaflucs_screen(K, M, Cnsq, dz, kz, r_filter)
%------------------------------------------------------

  iunit = i;
 
  dk = K(1,2);
  N = length(K);

  a1 = rand(N,N);
  a2 = rand(N,N);
  
  ksi0_flat = sqrt(3/2)/dk;
  ksi0_gauss = 1/sqrt(2)/dk;

  %ksi = ksi0_flat * ((2*a1-1) + iunit*(2*a2-1));
  ksi = ksi0_gauss * sqrt(-2*log(a1)) .* exp( iunit * 2*pi * a2) ;

  ksi = bereal(ksi);
  
  Phi = 0.033 * Cnsq * K.^(-11/3);
  Phi(1,1)=0;

  Sk = - ksi .* sqrt(Phi * 2*pi * dz) * kz;

  if (r_filter > 0)
    k_filter = 2*pi/r_filter;
    ind = find(K > k_filter);
    Sk(ind) = 0;
  end

  %Sk = ksi;

  %-- "continuous" Fourier transport back
  S = real(ifftn( Sk.* M*(N*dk)^2 ));

end

%------------------------------------------------------


function c = bereal(a)
% takes matrix of harmonics of complex matrix "a"
% and returns matrix of harmonics of a real matrix "c"
% assumes "a" is square of even size

  N = length(a);

  b = fftshift(a);

  b(N/2+2:N,2:N)   = flipud(fliplr(  conj(b(2:N/2,2:N))   ));

  b(N/2+1,N/2+2:N) =        fliplr( conj(b(N/2+1,2:N/2))  );
  b(1,N/2+2:N)     =        fliplr( conj(b(1,2:N/2))  );
  b(N/2+2:N,1)     =        flipud( conj(b(2:N/2,1))  );

  b( N/2+1, N/2+1 ) = real(b( N/2+1, N/2+1 ));
  b( 1,     N/2+1 ) = real(b( 1,     N/2+1 ));
  b( N/2+1, 1 )     = real(b( N/2+1, 1 ));
  b( 1, 1 )         = real(b( 1, 1 ));

  c = fftshift(b);

end

%------------------------------------------------------
