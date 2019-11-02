function guide_pearson_model

iunit=i;

 
  A0 = 60; t=0:0.00001:0.01;
  a0 = A0*[0.025, 0.1, 0.2, 0.26, 0.27, 0.4];  % initial amplitudes a0=sqrt(n)
    a0 = A0*0.3;

    for i = 1
 
     if (i==1) 
          hold off;   
     else  
          hold on;  
     end

     f0 = [0, (a0(i)).^2];

     f = lsode('fdot', f0, t);

     x = f(:,1); y = f(:,2);

     plot(x/pi,sqrt(y)/A0, 'b'); axis([0,2]);

     %plot(t,sqrt(y)/A0, 'b')
     %plot(t,x/pi)

  end

%return

  fid=fopen('model_A60_a0300.txt', 'wt');
  fprintf(fid, '# 1.t  2.amplitude/A0  3.theta/pi\n');
  for i=1:length(t)
	  fprintf(fid, '%16.6e  %16.6e  %16.6e\n', t(i), sqrt(y(i))/A0, x(i)/pi);
  end
  fclose(fid);


end

%---------------------

function f = fdot(x,t)

  %A=19.75;
  A=60;
  kk=1;

  N = A^2;

  f = [0 0]'; %'

  theta=x(1);
  n=x(2);
  
  f(1) = 2*kk  + 2*(N - 3*n)  + 2*(N - 4*n)*cos(theta);
  f(2) = 2*n*(N - 2*n) * sin(theta);

% f(1) = 2*kk  + 2*(N - 4*n)*(1 + cos(theta)) + 2*n;
% f(2) = 2*n*(N - 2*n) * sin(theta);
 
end

%---------------------


