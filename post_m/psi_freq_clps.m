function psi_freq_clps(fname_in, fname_out, N, L, kmax)
%
%       psi_freq_clps(fname_in, fname_out, N, L, kmax)
%
%   Compute correlation function of phase. Corelation function 
%   is either angle ageraged, or computed in x,y, xy, or yx directions.
%
%   Input:
%        fname_in       binary input file with psi
%        fname_out      text output file
%        N              size of simulation in points
%        L              domain size 
%        kmax           extend of integration in points
%
%----------------------------------------------------------------

  N0    = 11.7;

  f  =  read_psi(fname_in, N);
  dx = L/N;
  f0 = conj(f);

  Fm2 = zeros(N,N);
  Fm1 = zeros(N,N);
  F0 =  zeros(N,N);
  F1 =  zeros(N,N);
  F2 =  zeros(N,N);


  kkmax = kmax*kmax;

  for i=-kmax:kmax
     for j=-kmax:kmax

        kk = i*i + j*j;
        k  = sqrt(kk);

        if (kk > kkmax) 
	  continue;
        end
      
        dF = f0.*circshift(f,[-j,-i]); 

     if (k>0) 
        Fm2 = F2 + dF/kk; 
        Fm1 = F2 + dF/k; 
        F0 = F0 + dF; 
        F1 = F2 + dF*k; 
        F2 = F2 + dF*kk; 
     end

     end
  end

  Fm2 = Fm2;
  Fm1 = Fm1*dx;
  F0 = F0*dx^2;
  F1 = F1*dx^3;
  F2 = F2*dx^4;

  p = 0.5*real(F0);               % number of particles
  r = sqrt(real(F2)./real(F0));   % bump radius
 
  p1    = circshift(p,[0,1]); 
  p2    = circshift(p,[0,-1]); 
  p3    = circshift(p,[1,0]); 
  p4    = circshift(p,[-1,0]);

  ind   = find(p>p1 & p>p2 & p>p3 & p>p4 & p > N0); 
  [i,j] = find(p>p1 & p>p2 & p>p3 & p>p4 & p > N0); 

  p     = p(ind);
  r     = r(ind);
  h0    = abs(f(ind));

  h  = sqrt( p*2/pi )./r; 
  tc = 0.16*r.*r.*(p/N0 - 1).**(-0.62);
 
  freq = sum(1./tc);


%---------------------------

  Qm2 = real(Fm2(ind));
  Qm1 = real(Fm1(ind));
  Q0  = real(F0(ind));
  Q1  = real(F1(ind));
  Q2  = real(F2(ind));



    fid = fopen(fname_out,'wt');

    for k=1:length(p);
        fprintf(fid, '%8.6f  %12.6f  %12.6e  %12.6e  %12.6e  %12.6e\n', ...
		h0(k), Qm2(k), Qm1(k),  Q0(k),  Q1(k), Q2(k) );
    end

    fclose(fid);


  return;


  %-- writing output --


    fid = fopen(fname_out,'wt');
 
    fprintf(fid,'%% Expected frequency of collapses:  %12.6e\n', freq);
    fprintf(fid,'%% computed by \"psi_freq_clps.m\" from \"%s\".\n', fname_in);
    fprintf(fid,'%%\n%% 1.i  2.j  3.N/N0  4.r  5.h  6.h0  7.tc \n\n');

    for k=1:length(p);
        fprintf(fid, '%6d  %6d   %8.6f  %8.6f  %12.6e  %8.4f  %8.4f  %12.6e \n', ...
		i(k), j(k), p(k)/N0,  q(k),  r(k), h(k), h0(k), tc(k) );
    end

    fclose(fid);

  %-------------------------------------------------------------------



end

