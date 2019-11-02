function guide_oscil

  fbase='a064_cnd'; N = 256;

  NN=N*N;

  j1=N/2; j2=N/2+2;

  dt=0.001;

  %fid = fopen('oscil_nk.txt', "wt");
  fid  = fopen('oscil_psik.txt', "wt");

for k=1:1000
 
    num = num2str(k,'%04d');

    %-- amplitudes and phases --

    f = read_psi([fbase, '.psi.', num], N);
    f1=fft2(f)/NN; 
    f1=fftshift(f1);
    %f1 = f1.*conj(f1);

    fc = f1(j1:j2, j1:j2);
    fk(:,:,k) = fc;
    %f1(j1+1,j1+1) = 0;
    %nover = sum(sum(f1));


    fprintf(fid, "%8.4f", k*dt);
    fprintf(fid, "  %12.6e", abs(fc));
    %fprintf(fid, "  %12.6e", nover);
    fprintf(fid, "\n");


    %f2=flipud(fliplr(f1));
    %f2=circshift(f2,[1,1]);
    %sk = f1.*f2;

    %n0 = nk(N/2+1, N/2+1);
    %n5 = nk(j1:j2, j1:j2); n5(3,3)=0; n5=sum(sum(n5))/24;

    %imagesc(log(nk))
 
  end


fclose(fid);

end


%---------------------------------------------------------

function one_image(a, amin, amax, fname)

  a = max(a,amin);
  a = min(a,amax); 

  map = colormap(jet(256));

  sc   =  uint8(255*(a-amin)/(amax-amin)) + 1;
  imwrite(sc, map, fname);

end


%---------------------------------------------------------

function s = read_fft(fname, N)

   fid=fopen(fname, "rb");  
   f = fread(fid, N*N, 'double');  
   fclose(fid);  
   s = fftshift(reshape(f,N,N));

end


%---------------------------------------------------------
