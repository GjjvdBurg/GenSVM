#function test_kernel_dot_sparse(kerneltype)

  kerneltype = 'rbf'; 
 
  X = [1 2 0 0 0 0;
       0 3 0 4 0 0;
       0 0 5 6 7 0;
       0 0 0 0 0 8];
  [n, m] = size(X);
  spZ = sparse([ones(n, 1), X]);
  
  gamma = 0.05;
  degree = 1.7;
  const = 0.75;
  
  K = zeros(n, n);
  if strcmp(kerneltype, 'rbf')
    for i=1:n
      for j=1:n
        K(i, j) = exp(-gamma * sum((spZ(i, 2:end) - spZ(j, 2:end)).^2));
      end
    end
  end

#end