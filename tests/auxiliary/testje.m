clear;

rand('state', 123456);

n = 8;
m = 3;
K = 3;

y = [   2   1   3   2   3   3   1   2]';

U = SimplexGen(K);

UU = zeros(n, K-1, K);
 
for jj=1:K
  UU(:, :, jj) = U(y, :) - U(jj*ones(n, 1), :);
end

VV = zeros(n, K-1, K);
for i=1:n
  for j=1:K-1
    for k=1:K
      VV(i, j, k) = U(y(i), j) - U(k, j);
    end
  end
end

Z = [ones(n, 1), -1 + 2 * rand(n, m)];

V = -1 + 2 * rand(m+1, K-1);

ZV = Z*V;

Q = zeros(n, K);
for i=1:n
  for j=1:K
    Q(i, j) = ZV(i, :) * (U(y(i), :) - U(j, :))';
  end
end

% calculate loss
kappa = 0.5;
p = 1.5;
%rho = ones(n, 1);
rho = zeros(n, 1);
for i=1:K
  nk = sum(y == i);
  rho(y==i) = (n/(K*nk));
end
lambda = 0.123;

H = zeros(n, K);
for i=1:n
  for j=1:K
    q = Q(i, j);
    if (q <= -kappa)
      H(i, j) = (1 - q - (kappa + 1)/2.0);
    elseif (q <= 1)
      H(i, j) = (1/(2*kappa + 2)) * (1 - q)^2;
    else
      H(i, j) = 0;
    end
  end
end

R = zeros(n, K);
I = eye(K);
for i=1:n
  R(i, :) = I(y(i, :), :);
end
R = ~logical(R);

J = eye(m+1);
J(1, 1) = 0;

L = sum((H.^p).*R, 2).^(1/p);
L = 1/n * sum(rho.*L) + lambda * trace(V'*J*V);

% DON"T REMOVE YET!!