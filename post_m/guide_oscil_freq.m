
function guide_oscil_freq

  load 't60.dat';

  f = t60(:,7);
  tmax = 10;
  nmax = 10000;
  n1=150; n2=200;   % mode interval to search for max
  nstop = 500;      % max mode for plotting

  outfile = 't60_freq.dat';

  %----------------------

  fk=fft(f)/nmax; 
  ff=fk.*conj(fk);
  f2=flipud(ff); 
  f2=f2(1:end/2);


  [q,m]=max(f2(n1:n2)); m=m+n1-1,  w=2*pi*m/tmax

  plot((1:nstop),f2(1:nstop)); 


  fid = fopen(outfile, "wt");

  fprintf(fid, "%%1.mode  2.omega  3.fourier\n");

  for i=1:nstop
    fprintf(fid, "%6d  %12.6e %12.6e\n", i-1, 2*pi*(i-1)/tmax, f2(i));
  end

  fclose(fid); 
   

end


%---------------------------------------------------------
