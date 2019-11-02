function guide2d
% solves NLS equation 
%
%   i dPsi/dt + del^2(Psi) - |Psi|^2*Psi = 0 
%
% using 2nd order split-step method

%-- Input parameters --  

  fbase = 'a20_kx10_ky40';             % output file

  A0    =  20;              % \
  A1    =  0.1;              % -  amplitude of first harmonics
  A2    =  0.1;              % /

  L     =   2*pi;           % length of the box 
  N     =    256;           % discretization points 

  dt    =   1.e-4;           % timestep
  nmax  =   1.e+4;           % number of timesteps
  nplot =    1.e3;           % frequency of plotting
  nsave =     100;           % frequency of saving

  iunit = i;
  NN    = N*N;

%-- Initializaion --

  printf("\n Time step must be much smaller that period of phase rotation:\n");
  printf(" dt = %e,   T = %e.\n\n",  dt, 2*pi/A0^2);

  fU = zeros(N,N);

  %fU(1,2) = A1;  fU(2,1)   = A1;  fU(1,end) = A1;  fU(end,1)   = A1;
  %fU(2,2) = A2;  fU(2,end) = A2;  fU(end,2) = A2;  fU(end,end) = A2;
  %fU      =  fU*iunit/2;

  fU(40,10) = iunit;
  fU(1,1) = A0;

  U = ifftn(fU*NN);

  disp("time   N_spectral   N_real")
  disp([0,  sum(sum( fU.*conj(fU) )),  sum(sum( U.*conj(U) )) / N/N])


%-- save IC as binary file --

   p  = zeros(2,N,N);
   p(1,:,:) = real(U);
   p(2,:,:) = imag(U);

   p=reshape(p,N*N*2,1);

   fid = fopen([fbase,'.psi.0000'], 'wb');
   fwrite(fid, p, 'double');
   fclose(fid);

%-- precompute exponents for linear step --

  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);
  [kx,ky] = meshgrid(k,k);

  DD = -(kx.^2 + ky.^2)*(2*pi/L)^2;
 
  ff = exp(-(DD/(0.24*N)^2).^64);

  %ind=find(ff==1); ff(ind)=2;

  %plot(k,ff(1,:))
  %imagesc(fftshift(ff)); axis('equal')
  %max(max(ff))


  w=1/(2 - 2^(1/3));
 
  expLin1 = exp(iunit*DD*dt*w);
  expLin2 = exp(iunit*DD*dt*(1-2*w));
  expLin3 = expLin1;


%-- non-linear timesteps --

  dt0= 0.5*w*dt;
  dt1= 0.5*(1-w)*dt;
  dt2= 0.5*(1-w)*dt;
  dt3= 0.5*w*dt;

  izero = N/4:3*N/4+1;  


%-- Main loop - time advancement --

  fid = fopen ([fbase, '.txt'], 'wt');

  fprintf(fid, '#1.t  2.Nwaves  3.(0,0)   '); 
  fprintf(fid, '4.(1,0) 5.(-1,0)   6.(0,1) 7.(0,-1)  '); 
  fprintf(fid, '8.(1,1) 9.(-1,-1)  10.(1,-1) 11.(-1,1)\n\n'); 

  for n=0:nmax

    %------------------------

    %--advance non-linear part 0

    UU = U.*conj(U);
    U  = U .* exp(-iunit*dt0*UU);


    %-- advance linear part 1

    fU = fftn(U/NN);            % FT to k-space
    fU = fU.*ff;
    fU = fU .* expLin1;        % time advancement
    fU(izero, :) = 0;          % dealiasing
    fU(:, izero) = 0;          % dealiasing
    U  =  ifftn(fU*NN);         % FT to r-space


    %-- advance non-linear part 1

    UU = U.*conj(U);
    U  = U .* exp(-iunit*dt1*UU);

    %-- advance linear part 2

    fU = fftn(U/NN);            % FT to k-space
    fU = fU .* expLin2;        % time advancement
    fU(izero, :) = 0;          % dealiasing
    fU(:, izero) = 0;          % dealiasing
    U  =  ifftn(fU*NN);         % FT to r-space
 
    %-- advance non-linear part 2

    UU = U.*conj(U);
    U  = U .* exp(-iunit*dt2*UU);


    %-- advance linear part 3

    fU = fftn(U/NN);            % FT to k-space
    fU = fU .* expLin3;        % time advancement
    fU(izero, :) = 0;          % dealiasing
    fU(:, izero) = 0;          % dealiasing
    U  =  ifftn(fU*NN);         % FT to r-space
 
    %--advance non-linear part 3

    UU = U.*conj(U);
    U  = U .* exp(-iunit*dt3*UU);

    %------------------------
    %-- plot and save

    if (rem(n,nsave) == 0) 
 
      fU = fftn(U/NN);
      a  = abs(fU);

      nWaves = sqrt(sum(sum(fU.*conj(fU))));
      fprintf(fid, '  %10.4e   %12.6e   %12.6e',  n*dt,  nWaves, a(1,1) );
      fprintf(fid, '   %10.4e %10.4e   %10.4e %10.4e', ...
	             a(1,2), a(1,end), a(2,1), a(end,1) );
      fprintf(fid, '   %10.4e %10.4e   %10.4e %10.4e\n', ...
	             a(2,2), a(end,end),  a(2,end), a(end,2) );

    endif

    if (rem(n,nplot) == 0) 
      imagesc(flipud(fftshift(log(abs(fU))))); 
      axis('equal'); axis('tic','off');

      disp([n*dt, sum(sum( fU.*conj(fU) )), sum(sum( U.*conj(U) )) / N/N])

      fidbin = fopen([fbase,'.fk.', num2str(n/nplot, '%04d')], 'wb');
      fwrite(fidbin, abs(fU), 'double');
      fclose(fidbin);

      pause(0.2);
    endif


    %------------------------
 
  end

  fclose(fid);

  

end


%-------------------------------------------------

