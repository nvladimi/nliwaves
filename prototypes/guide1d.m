function guide1d
% solves NLS equation 
%
%   i dPsi/dt + del^2(Psi) - |Psi|^2*Psi = 0 
%
% using 2nd order split-step method

%-- Input parameters --  

  %fileout = "tmp.txt";     % output file
  fileout = "o4_dt1e-4_n0256_short.txt";     % output file

  A0    = 100;              % \
  A1    =   0;              % -  amplitude of first harmonics
  A2    =   1;              % /

  L     =   2*pi;           % length of the box 
  N     =    256;           % discretization points 

  dt    =   1.e-4;           % timestep
  nmax  =   1.0e3;           % number of timesteps
  nplot =   1.0e3;           % frequency of plotting
  nsave =       1;           % frequency of saving

  iunit = i;


%-- Initializaion --

  printf("Time step must be much smaller that period of phase rotation:\n");
  printf("dt = %e,   T = %e.\n\n",  dt, 2*pi/A0/A0);

  fU = zeros(1,N);

  fU(1) = A0;
  fU(2) = iunit*A1/2;  fU(end)   = iunit*A1/2;
  fU(3) = iunit*A2/2;  fU(end-1) = iunit*A2/2;

  fU    = fU*N;


  U = ifft(fU);

  %plot(1:N, abs(fU));


%-- precompute exponent for linear step

  k = 0:N-1; 
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);

  DD = -(k*2*pi/L).^2;

  %semilogy(k, abs(fU/N))

  w=1/(2 - 2^(1/3));
 
  expLin1 = exp(iunit*DD*dt*w);
  expLin2 = exp(iunit*DD*dt*(1-2*w));
  expLin3 = expLin1;


%-- non-linear timesteps

  dt1= 0.5*(1-w)*dt;
  dt2= dt1;
  dt3= w*dt;


%-- Main loop - time advancement --

  ind=1:2:13; 
  format short e


  fid = fopen (fileout, 'wt');

  fprintf(fid, '#1.t 2.sqrt(N)  3.k=0,2,4,6,8,10,12\n\n'); 

  for n=0:nmax


    %------------------------

    % advance linear part 1

    fU = fft(U);               % FT to k-space
    fU = fU .* expLin1;        % time advancement
    fU(N/4:3*N/4+1) = 0;       % dealiasing
    U  =  ifft(fU);            % FT to r-space
 
    % advance non-linear part 1

    UU = U.*conj(U);
    U  = U .* exp(-iunit*dt1*UU);


    % advance linear part 2

    fU = fft(U);               % FT to k-space
    fU = fU .* expLin2;        % time advancement
    fU(N/4:3*N/4+1) = 0;       % dealiasing
    U  =  ifft(fU);            % FT to r-space
 
    % advance non-linear part 2

    UU = U.*conj(U);
    U  = U .* exp(-iunit*dt2*UU);


    % advance linear part 3

    fU = fft(U);               % FT to k-space
    fU = fU .* expLin3;        % time advancement
    fU(N/4:3*N/4+1) = 0;       % dealiasing
    U  =  ifft(fU);            % FT to r-space
 
    % advance non-linear part 1

    UU = U.*conj(U);
    U  = U .* exp(-iunit*dt3*UU);


    %------------------------

    % plot and save

    if (rem(n,nsave) == 0)

      nWaves = sqrt(sum(fU.*conj(fU)))/N;
      fprintf(fid, ' %10.4e   %10.4e   ', n*dt, nWaves, abs(fU(ind))/N );
      fprintf(fid, ' %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n', ...
	            abs(fU(ind))/N );
   endif

    if (rem(n,nplot) == 0) 


      %disp(abs(fU(ind)/N))
	 %semilogy(k, abs(fU/N) );  axis([-N/4, N/4, 1e-12, 1.e2] )
	   bar(log(abs(fU/N))); axis([0, N, -30, 5] ); set (gca, "xtick", [0:32:N]); 

      pause(0.1)
    endif


    %------------------------
 
  end

  fclose(fid)

end


%-------------------------------------------------

