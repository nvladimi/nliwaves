function twomode_post(fbase, fnums)
 
    fname = [fbase, '.param'];
  
    load(fname);   %('fbase', 'fnum', 'epsilon', 'gamma', 'theta', 'n'); 
  
    N = n(3)*n(4);

    A = zeros(N*length(fnums), 5);

    i=0;

    for fnum = fnums 

        fname    = [fbase, '.',  num2str(fnum, '%04d')];

        fid = fopen([fname, '.a1a2'], 'rb');

        a = fread(fid, N*5, 'double');
        a = reshape(a, [N, 5]);

        fclose(fid);

        A(i*N+1: (i+1)*N, :) = a;
  
        i = i+1;

    end

%-- debugging -- 


    A = reshape(A, [N*length(fnums), 5]);
    t = A(:,1);
    A = A(:, 2:5);
    

    [N1, N2, Hnl] = hamiltonian(A);


   figure(1);
   plot(t, A(:,1), '-r', t, A(:,2), '--r', t, A(:,3), '-b',  t, A(:,4), '--b' );
axis([0, 10*pi, -20,20]); grid("on");
   set(gca, "fontsize", 20);

   figure(2);
plot(t, N1, '-g', t, N2, '-m',   t, N1+2*N2, '-b',  t, Hnl, '-r');
axis([0, 2*pi, -200,500]); grid("on");
   set(gca, "fontsize", 20)


end

%---------------------


%---------------------

 function [N1, N2, Hnl] = hamiltonian(x);

  global epsilon;

  N1  = x(:,1).*x(:,1) + x(:,2).*x(:,2);
  N2  = x(:,3).*x(:,3) + x(:,4).*x(:,4) ;
  Hnl  = epsilon * 2 * ( x(:,1).*x(:,1).*x(:,3) - x(:,2).*x(:,2).*x(:,3) + 2*x(:,1).*x(:,2).*x(:,4) ) ;

 
end

%---------------------
