
function clps_pdfmax(fname_in, fname_out)
%
%   clps_pdfmax(fname_in, fname_out)
%
%   Collecting PDF of collapse maxima from grep'ed headers
%
%   Input: name of file with data grep'ed from  *.post.?? or *.clps.?? 
%         "post" files are better because they exclude abrupted collapses.
%
%   grep collapse *.post* | awk '{printf("%12.6e %12.6e\n", $4, $6)}' > all.hmax 


  %----------------------------------------------

  %-- collect only collpses above 6

  allh = importdata(fname_in);

  h = allh(:,2);  i=h>6;  h = allh(i,:);


  %-- remove dublicate entries

  [t,h] = remove_doubles(h);


  %-- just in case we want to list all collpases 
  % 
  % fid = fopen('tmp.dat','wt');
  % for k=1:length(h)
  %     fprintf(fid, '%4d  %10.4f  %10.4f\n', k, t(k), h(k));
  % end
  % fclose(fid);


  %-- collect distributions

  b1  = 0.5:1:1000;       % bins size of 1
  d1  = hist(h, b1);

  b10 = 5:10:1000;        % bins size of 10
  d10 = hist(h, b10);


  %-- print distribution to file

  fid = fopen(fname_out,'wt');
  
  fprintf(fid, '%% Created by \"clps_pdfmax.m\" from greped headers.\n');
  fprintf(fid, '%% Dataset: %d collapses with hmax > 6.\n', length(h));
  fprintf(fid, '%% First output set is for bins of size 10,\n');
  fprintf(fid, '%% second set is for bins of size 1.\n\n');
  fprintf(fid, '%% 1.h_bin_centers  2.collapses_in_bin  3.bins_integrated\n\n');

  p=0;
  for k=length(d10):-1:2
     p =p + d10(k);
     fprintf(fid, '%8.2f  %10d  %10d\n', b10(k), d10(k), p);
  end

  fprintf(fid, '\n\n');
 
  p=0;
  for k=length(d1):-1:7
     p =p + d1(k);
     fprintf(fid, '%8.2f  %10d  %10d\n', b1(k), d1(k), p);
  end

  fclose(fid);

end

%-----------------------------------------------------

function [t,h] = remove_doubles(h2)

  [h,i]   = sort(h2(:,1));
  h2(:,1) = h;
  h2(:,2) = h2(i,2);

  for i=2:length(h2)
    if ( h2(i,:)== h2(i-1,:))
      h2(i-1,2) = NaN;
    end
  end

  h=h2(:,2);
  i=isfinite(h);
  h=h(i);
  t=h2(i,1);

end
