function clps_betatau_filter(fname_in, fname_out, threshold, betamax)
%
%   clps_betatau_filter(fname_in, fname_out, threshold)
%
%   Replaces beta_t with beta_tau; skips data at the refinement jumps.
%


  %-- read data in --

  data = importdata(fname_in);

  %-- replace 7th column (beta_t) with beta_tau

  data(:,7)=data(:,7).*data(:,4).^2;

  b=data(:,7);
  a=data(:,6); 

  db = (circshift(b,[-1,0])-b).*(b-circshift(b,[1,0]));
  db = db./b.^3 ;

  ind = find ( abs(db) < threshold | a > betamax );
  filtered_data = data(ind,:);


  %-- output --

  fid=fopen(fname_out, 'wt');

  fprintf(fid, ... 
    '%% created by \"clps_betatau_filter.m\" from \"%s\"\n\n', fname_in);
  fprintf(fid, ... 
    '%%1.t  2.h  3.tmax-t  4.L  5.|psi|_rr  6.beta  7.beta_tau  ');
  fprintf(fid, ... 
    '8.phase  9.phase_rr  10.V  11.ss_time  12.ss_a\n');
 


  for i=1:length(filtered_data) 
      fprintf(fid, ...
	' %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e  %22.15e\n',...
      filtered_data(i,:));
  end

  fclose(fid);


end



end
