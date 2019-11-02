function contour_lenght

fname_in = 'contour_m4_a1e2_0005_f1e2_0005.0008';
fname_out  = 'lenght_m4_a1e2_0005_f1e2_0005.0008';


%c = importdata(fname_in);

load c.txt

l = circshift(c,[-1,;]) - c;
l = l(1:end-1,:);

l = sqrt( l(:,1).*l(:,1) + l(:,2).*l(:,2) ); 
sum(l)

for s=2:length(l)
	l(s) = l(s) + l(s-1)
end

r=c;  r0=c(1,:);

r(:,1) = r(:,1) - r0(1);
r(:,2) = r(:,2) - r0(2);

r = sqrt( r(:,1).*r(:,1) + r(:,2).*r(:,2) );
r=r(2:end);




  % fid=fopen(fname_out, 'at');
  % fprintf(fid,  '%6d  %8.4f  %13.6e\n', M, dx, L);
  % fclose(fid);

end

