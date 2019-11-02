function nls_ic

  %fname_out='g05h10r113m10.psi.0000';
  %N=512;  L=12.8; h0=1.0;  r0=0.61;  nc=9; Mcorner = 1.0;
  %N=8; L=0.8; nc=1;

   fname_out='test.psi.0000';
   N=512;  L=25.6; h0=0.5;  r0=1.13;  nc=10; Mcorner = 1.0;



  rand(3);

  psi = zeros(N,N);

  %--- grid of gaussians with random phases

  for ic=0:nc-1

    x0 = (ic+0.5)*L/nc;
    x = (0:N-1)*L/N + x0;

    ind=find(x<-L/2); 
    x(ind)=x(ind)+L;

    ind=find(x>L/2); 
    x(ind)=x(ind)-L;

    for jc=0:nc-1

      y0 = (jc+0.5)*L/nc;

      y = (0:N-1)*L/N + y0;

      ind=find(y<-L/2); 
      y(ind)=y(ind)+L;

      ind=find(y>L/2); 
      y(ind)=y(ind)-L;

      [X, Y] = meshgrid (x, y);
 
      a0 = h0*exp(rand*2*pi*i);
      h = a0*exp( - (X.^2 + Y.^2)/r0^2 );

      psi = psi + h;

    end

  end

  
  %--- domain-size modulation by cos


    x0 = 0.5*L;
    x = (0:N-1)*L/N + x0;

    ind=find(x<-L/2); 
    x(ind)=x(ind)+L;

    ind=find(x>L/2); 
    x(ind)=x(ind)-L;


    y0 = 0.5*L;
    y = (0:N-1)*L/N + y0;

    ind=find(y<-L/2); 
    y(ind)=y(ind)+L;

    ind=find(y>L/2); 
    y(ind)=y(ind)-L;

    [X, Y] = meshgrid (x, y);
 
    M = cos(2*pi*X/L)+cos(2*pi*Y/L);
    M = (M+2)/4*(1-Mcorner) + Mcorner;
 
    psi = psi.*M;


  %-- trim specrum to avoid refinement

    p = fftn(psi);
    ind=N/8:N-N/8;
    p(ind,:) = 0;
    p(:,ind) = 0;
    psi = ifftn(p);

  %-- save IC to file

    N_particles = sum(sum(psi.*conj(psi)))*(L*L)/(N*N)
    per_collapse = N_particles/(nc*nc)
    avg_psi = sum(sum(abs(psi)))/(N*N)

    imagesc(abs(psi))


    %p=fftn(psi)/N/N;    
    %pp = p.*conj(p);
    %plot(1:N, pp(1,:))
    %imagesc(log(pp));

    return

    p  = zeros(2,N,N);
    p(1,:,:) = real(psi);
    p(2,:,:) = imag(psi);

    p=reshape(p,N*N*2,1);

    f = fopen(fname_out, 'wb');
    fwrite(f, p, 'double');
    fclose(f);

end

