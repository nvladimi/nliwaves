function rebin_pdf(fname_in, fname_out, rebin)
%
%   rebin_pdf(fname_in, fname_out, rebin)
%
%   resample PDFs into larger bins
%
%   Input:
%       fname_in    test input file
%       fname_out   text output file
%       rebin       groups of bins combined 


  p = importdata(fname_in);

%load(fname_in); p=loc_max_hist;


  nbins = length(p)/rebin;

  %-- output --

  fid = fopen(fname_out,'wt');
  fprintf(fid,'%% Probability density function (PDF) of Psi \n');
  fprintf(fid,'%% computed by \"rebin_pdf.m\" from \"%s\".\n', fname_in);
  fprintf(fid,'%%\n%% 1.bin  2.psi  3.pdf(abs)  4.pdf(Re)  5.pdf(Im)\n\n');

  for k=1:nbins

    i1 = (k-1)*rebin+1;
    i2 = k*rebin;
    pdf = sum(p(i1:i2,:))/rebin;

%    fprintf(fid,'%6d %8.4f  %16.8e  %16.8e  %16.8e\n', ... 
%                k, pdf(2), pdf(3), pdf(4), pdf(5));

    fprintf(fid,'%8.4f  %16.8e  %16.8e  %16.8e\n', ...
                pdf(2), pdf(3), pdf(4), pdf(5));


%    fprintf(fid,'%8.4f  %16.8e\n', pdf(1), pdf(2));

  end

  fclose(fid);

end
