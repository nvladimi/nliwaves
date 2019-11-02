
n = 1024*4

fbase_out ='driving/all_4k.map'
fbase_in  = 'driving/m4_a1e2_f1e2.map'

%------------------------

s=[];

D = contour_map_avg([fbase_in,'.0001'], [fbase_out,'.0001'], n, [0:57]);
s = [s, D(end,:)];

D = contour_map_avg([fbase_in,'.0002'], [fbase_out,'.0002'], n, [0:52]);
s = [s, D(end,:)];

D = contour_map_avg([fbase_in,'.0003'], [fbase_out,'.0003'], n, [0:48]);
s = [s, D(end,:)];

D = contour_map_avg([fbase_in,'.0004'], [fbase_out,'.0004'], n, [0:63]);
s = [s, D(end,:)];

D = contour_map_avg([fbase_in,'.0005'], [fbase_out,'.0005'], n, [0:57]);
s = [s, D(end,:)];

D = contour_map_avg([fbase_in,'.0006'], [fbase_out,'.0006'], n, [0:47]);
s = [s, D(end,:)];



reavg_any(fbase_out,  [fbase_out,'.all'], [1:6] ) 


avg = sqrt(sum(sum(s.*s))/length(s));
[y,x] = hist(s, [-2000:20:2000], 1/20. );

gx=-5:0.1:5; gy=1/sqrt(2*pi)*exp(-gx.*gx/2);  % gauss
semilogy(x/avg,y*avg, '*', gx, gy);


%----------

fid = fopen([fbase_out, '_pdf.txt'], 'wt');

fprintf(fid, '%% created by \"driving_pdf.m\" from \"%s.*\"\n', fbase_in);
fprintf(fid, '%% Standart deviation:  %16.8e\n', avg);
fprintf(fid, '%% 1.ksi/avg  2.pdf*avg\n\n');

for i=1:length(x)
  fprintf(fid,' %16.8e  %16.8e\n', x(i)/avg,  y(i)*avg);
end

fclose(fid);




