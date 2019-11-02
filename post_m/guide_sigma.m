function f5=guide_sigma

  fbase='a160e-2s'; N = 512;
  %fbase='a064s'; N = 256;
  %fbase='tmp'; N = 256;

  NN=N*N;

  nkmin  = -36;
  nkmax  =  6;


  fftot = zeros(256,256);

  i1=(1:N);
  i2=[1,(N:-1:2)];


  %j1=N/2; j1=N/2+2;
  j1=N/2-1; j2=N/2+3;


  q=[0,0,0,0,0,0, 0];

for k=10:11 %:1:1061
 
    num = num2str(k,'%04d');

    %-- amplitudes and phases --
[fbase, '.psi.', num]

    f = read_psi([fbase, '.psi.', num], N);
    f1=fft2(f)/NN; 
    f1=fftshift(f1);

    f2=flipud(fliplr(f1));
    f2=circshift(f2,[1,1]);

    ff = f1.*f2;
    nk = f1.*conj(f1);

    n  = sum(sum(nk));
    n0 = nk(N/2+1, N/2+1);
    f0 = ff(N/2+1, N/2+1);

    phi=(arg(ff(j1:j2, j1:j2)) - arg(ff(N/2+1, N/2+1)) )/pi
    amp=abs(ff(j1:j2, j1:j2))

    n5 = nk(j1:j2, j1:j2);   % n5(3,3)=0; %n5=sum(sum(n5))/24;
    f5 = ff(j1:j2, j1:j2); f5(3,3)=0; f5=sum(sum(f5))/24;

    %t=0.001*k;
    %a=atan(524*t);

    %q=[q; [t, n,n0, n5,abs(f5),  arg(f0), arg(f5)]];

 
    continue


    for i=1:N
      ff1(i) = abs(ff(i,i));
      ff2(i) = abs(ff(i, i2(i)));
      nk1(i) = nk(i,i);
      nk2(i) = nk(i, i2(i));
    end

    i=-N/2:N/2-1;
    semilogy(i,ff1,'or',  i,ff2, '+r',  i,nk1,'ob', i,nk2,'xb');
    axis([-32,32, 1e-6, 1.e3]);
    %axis([-20,20, 0, 10]);

    one_image(abs(ff)-abs(nk), -1.e-2, 1.e-2, [fbase, '.diff.', num, '.png']);
    one_image(log(abs(ff)), nkmin, nkmax, [fbase, '.kkh.', num, '.png']);
    one_image(log(abs(nk)), nkmin, nkmax, [fbase, '.nk.', num, '.png']);
    one_image(arg(ff), -pi,   pi,    [fbase, '.kkp.', num, '.png']);
    %continue
 
  end

  q=q(2:end,:);


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
