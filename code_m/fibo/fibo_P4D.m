function fibo_P4D(fbasein, fbaseout, fstart, seeds, i0, nB, dB)
%
%  Input:
%
%   fbasein     - string with base name for input files
%   fbaseout    - string with base name for output files
%   fstart      - file number to start with, in each seed
%   fseeds      - array of seed numbers to analyze
%   i0          - index of the third mode in triplet to consider
%   dB          - width of the bins for n_i, in terms of <n_i>
%   nB          - size of 4D arrays, number of bins for modes i0, i0-1, i0-2, theta;
%

   global P4D;            % 4D probability array
   global P3D;            % 2D probability array
   global nBins;          % number of bins for modes i0, i0-1, i0-2, theta; 
   global dBins;          % width of the bins i0, i0-1, i0-2;
   
   nBins = nB;            % number of bins for modes i0, i0-1, i0-2, theta; 
   dBins = dB;            % width of the bins i0, i0-1, i0-2;

%-- read data and compute average --

   ntot = 0;

   aa_avg = zeros(1,3);

   for s = seeds

       fbase = [fbasein, '_s', num2str(s)];
       fbase = [fbase, '/', fbase];

       fnum = fstart;

       while 1

           fname    = [fbase, '.',  num2str(fnum, '%04d'), '.ak'];

           if exist(fname, "file")

              load([fbase, '.',  num2str(fnum, '%04d'), '.param']);
 	   
              fid = fopen(fname, 'rb');
              a = fread(fid, nsave*(2*m+1), 'double');
              a = reshape(a, [nsave, 2*m+1]);

              ind = (1:m)*2;
              a = a(:, ind)  + 1i*a(:, ind+1);

              a = a(:, i0-2:i0);                 % <-- a(nsave,3);
              aa = a.*conj(a);

              aa_avg = aa_avg + sum(aa,1);

              fnum = fnum+1;
              ntot = ntot + nsave;
 
           else % file does not exist
             break
           end

       end  % loop over files

   end % loop over seeds

   aa_avg = aa_avg/ntot;
		    
   fnameout = [fbaseout, '.param'];

   save( fnameout, 'i0', 'nBins', 'dB', 'ntot',  'aa_avg'); 

   dBins = dB*aa_avg;


%-- create probability array  --

   P4D = zeros(nBins(1), nBins(2),  nBins(3), nBins(4), 'uint32');
   P3D = zeros(nBins(2),  nBins(3), nBins(4), 'uint32');

%-- read data and assemble distribution --

   ntot = 0;

   for s = seeds

       fbase = [fbasein, '_s', num2str(s)];
       fbase = [fbase, '/', fbase];
 
       fnum = fstart;

       while 1

           fname    = [fbase, '.',  num2str(fnum, '%04d'), '.ak'];

           if exist(fname, "file")

              load([fbase, '.',  num2str(fnum, '%04d'), '.param']);
 	   
              fid = fopen(fname, 'rb');
              a = fread(fid, nsave*(2*m+1), 'double');
              a = reshape(a, [nsave, 2*m+1]);

              ind = (1:m)*2;
              a = a(:, ind)  + 1i*a(:, ind+1);

              a = a(:, i0-2:i0);                 % <-- a(nsave,3);

              aa = a.*conj(a);

              phi = angle(a(:,3)) - angle(a(:,2)) - angle(a(:,1));

              ind = find (phi >  pi);  phi(ind) = phi(ind) - 2*pi;
              ind = find (phi >  pi);  phi(ind) = phi(ind) - 2*pi;
              ind = find (phi < -pi);  phi(ind) = phi(ind) + 2*pi;
              ind = find (phi < -pi);  phi(ind) = phi(ind) + 2*pi;

              phi2 = angle(a(:,3)) - angle(a(:,2));

              ind = find (phi2 >  pi);  phi2(ind) = phi2(ind) - 2*pi;
              ind = find (phi2 < -pi);  phi2(ind) = phi2(ind) + 2*pi;
  
              add_counts(aa, phi, phi2);

              fnum = fnum+1;
              ntot = ntot + nsave;
 
           else % file does not exist
             break
           end

       end  % loop over files

   end % loop over seeds
  
   fnameout = [fbaseout, '.param'];
   save( fnameout, 'i0', 'nBins', 'dB', 'ntot',  'aa_avg'); 

   fid = fopen([fbaseout, '.P4D'], 'wb');
   fwrite(fid, P4D, 'uint32');
   fclose(fid);

   fid = fopen([fbaseout, '.P3D'], 'wb');
   fwrite(fid, P3D, 'uint32');
   fclose(fid);

return

end


%=======================================

  function add_counts(aa, phi, phi2)

    global P4D;
    global P3D;
    global nBins;
    global dBins;

    ntot = length(aa);

    nPhi = nBins(4);
    nPi  = nPhi/2;

    iN1 = floor( aa(:, 1) / dBins(1) ) + 1;
    iN2 = floor( aa(:, 2) / dBins(2) ) + 1;
    iN3 = floor( aa(:, 3) / dBins(3) ) + 1;

    iN1  = min(iN1, nBins(1));
    iN2  = min(iN2, nBins(2));
    iN3  = min(iN3, nBins(3));

    iPhi = round( phi./pi * nPi) + nPi + 1;
    ind = find(iPhi == nPhi+1);
    iPhi(ind) = 1;

    iPhi2 = round( phi2./pi * nPi) + nPi + 1;
    ind = find(iPhi2 == nPhi+1);
    iPhi2(ind) = 1;


    for i=1:ntot
	    P4D(iN1(i),iN2(i),iN3(i),iPhi(i)) = P4D(iN1(i), iN2(i), iN3(i), iPhi(i))  + uint32(1);
	    P3D(iN2(i),iN3(i),iPhi2(i)) = P3D(iN2(i), iN3(i), iPhi2(i))  + uint32(1);
    end
  
  end

%=======================================

