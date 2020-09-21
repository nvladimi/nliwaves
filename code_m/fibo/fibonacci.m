function fn = fibonacci(M)

    fn = zeros(1, M);
    fn(1) = 1;
    fn(2) = 1;
    for k=3:M
	fn(k) = fn(k-1)+fn(k-2);
    end
    
end

