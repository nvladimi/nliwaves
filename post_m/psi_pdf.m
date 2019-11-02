function psi_pdf(fname_in, fname_out, N, nbins, limits)
%
%   psi_pdf(fname_in, fname_out, N, nbins, limits)
%
%   compute PDF of abs(Psi), re(Psi), Im(Psi) 
%
%   Input:
%       fname_in    binary input file
%       fname_out   text output file
%       N           size of simulation
%       nbins       number of bins
%       limits      [psimin, psimax]

  f = read_psi(fname_in, N);

  u = real(f);
  v = imag(f);
  h = abs(f);

  u = reshape(u, N*N, 1);
  v = reshape(v, N*N, 1);
  h = reshape(h, N*N, 1);

  %-- find bin centers  

  dn = (limits(2)-limits(1))/nbins;
  bins = linspace(limits(1)+dn/2, limits(2)-dn/2, nbins);

  pU = hist(u, bins, 1);
  pV = hist(v, bins, 1);
  pH = hist(h, bins, 1);

  %-- output --

  fid = fopen(fname_out,'wt');
  fprintf(fid,'%% Probability density function (PDF) of Psi \n');
  fprintf(fid,'%% computed by \"psi_pdf.m\" from \"%s\".\n', fname_in);
  fprintf(fid,'%%\n%% 1.bin  2.psi  3.pdf(abs)  4.pdf(Re)  5.pdf(Im)\n\n');

  for k=1:length(bins)

    fprintf(fid,'%6d %8.4f  %16.8e  %16.8e  %16.8e\n', ...
	    k-1, bins(k), pH(k), pU(k), pV(k));

  end

  fclose(fid);

end
