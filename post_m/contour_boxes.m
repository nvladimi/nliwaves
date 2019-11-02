function contour_boxes(fnum, crange)

%fname_in='contours/m4_a1e2_f1e2.contraw.0006.0015';
%fname_out='boxes/m4_a1e2_f1e2.boxes.0006.0015';

fbase = 'm4_a1e2_f1e2';
N     = 4096;

%-------------------

boxes = 2.^(3:log2(N));
%boxes=512

for nc=crange 

   ext = [num2str(fnum,'%04d'), '.', num2str(nc,'%04d')];
   fname_in =  ['contours/', fbase, '.contraw.', ext];
   fname_out = ['boxes/',    fbase, '.boxes.',   ext];


   c_in = importdata(fname_in);

   %load c_in.txt

   for M=boxes

      dx = N/M;
      m = zeros(M,M);
      c = ceil(c_in/dx); 

      for (k=1:length(c)) 

         i=c(k,1);
         j=c(k,2);
         m(j,i) = 1;            % unweighted
         %m(j,i) = m(j,i) + 1;  % weigted by number of hits

      end

      %max(max(m))

      %imagesc(flipud(m)); axis('equal', 'off');
      L=sum(sum(m));

      fid=fopen(fname_out, 'at');
      fprintf(fid,  '%6d  %8.4f  %13.6e\n', M, dx, L);
      fclose(fid);

   end
end

