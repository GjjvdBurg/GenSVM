function newV = update_test
  rand('state', 2190382);

  more off;
  n = 8;
  m = 3;
  K = 3;
  
  p = 1.1;
  lambda = 0.123;
  kappa = 0.5;
  
  y = [   2   1   3   2   3   3   1   2]';
  
  
  U = SimplexGen(K);
  UU = zeros(n, K-1, K);
  for jj=1:K
    UU(:, :, jj) = U(y, :) - U(jj*ones(n, 1), :);
  end
  
  R = zeros(n, K);
  I = eye(K);
  for i=1:n
    R(i, :) = I(y(i, :), :);
  end
  R = ~logical(R);
  
  Z = [ones(n, 1), -1 + 2 * rand(n, m)];
  V = -1 + 2 * rand(m+1, K-1);
  rho = ones(n, 1);
  
  V
  
  newV = getUpdate(Z, y, rho, p, kappa, lambda, UU, R, V, false);
end

function set_vector(a)
  for i=1:numel(a)
    fprintf('A[%i] = %.16f;\n', i-1, a(i));
  end
end

function set_matrix(A)
  for i=1:size(A, 1)
    for j=1:size(A, 2)
      fprintf('matrix_set(A, %i, %i, %i, %.16f);\n', size(A, 2), i-1, j-1, A(i, j));
    end
  end
end

% From the original Matlab implementation
function [newV, oldL] = getUpdate(Z, y, rho, p, kappa, lambda, UU, R, V, debug)
  % Initialize constants
  [n, m] = size(Z);
  m = m - 1;
  [~, K] = size(R);
  
  % Calculate initial errors
  ZV = Z*V;
  q = zeros(n, K);
  for jj = 1:K
      q(:, jj) = sum(ZV.*UU(:,:,jj), 2);
  end
  
  % Calculate Huber hinge errors
  G1 = (q <= -kappa);
  G2 = (q <= 1) & (~G1);
  H = (1 - q - (kappa+1)/2).*G1 + (1/(2*kappa + 2))*((1 - q).^2).*G2;
  
  % Split objects in Category 1 and 2
  C = sum((H.*R)>0, 2)<=1;
  C1i = find(C>0); % for some reason this is faster than C==1
  C2i = find(C<1); % for some reason this is faster than C==0
  
  % We are going to do the calculations separately for each category (1,2)
  % so we separate all important variables as well.
  Q1 = q(C1i, :);
  Q2 = q(C2i, :);
  
  R1 = R(C1i, :);
  R2 = R(C2i, :);
  
  P1 = rho(C1i, :);
  P2 = rho(C2i, :);
  
  y1 = y(C1i, :);
  
  H1 = H(C1i, :);
  H2 = H(C2i, :);
  
  ZV1 = ZV(C1i, :);
  ZV2 = ZV(C2i, :);
  
  UU1 = UU(C1i, :, :);
  UU2 = UU(C2i, :, :);
  
  n1 = length(y1);
  n2 = n - n1;
  
  %% First create all matrices for the Case 1 objects.
  if debug
      % Check the first majorization of the Case 1 objects.
      TT1 = sum((H1.^p).*R1,2).^(1/p);
      TT2 = sum(H1.*R1,2);
      if abs(sum(TT1) - sum(TT2))/n1 > eps
          fprintf('\tFirst Case 1: diff = %15.16f\n', abs(sum(TT1) - sum(TT2))/n1);
      end
  end
  
  % First do the majorization for the Case 1 objects (p = 1 in Huber maj.)
  G1 = (Q1 <=-kappa);
  G2 = (Q1 <= 1) & (~G1);
  G3 = ~(G1|G2);
  
  % calculate dummy variables
  Phi = 1 - Q1 - (kappa + 1)/2;
  Psi = (1 - Q1)/sqrt(2*kappa + 2);
  iPhi = 1./Phi;
  
  a1 = 1/4*iPhi.*(G1 - G3) + (1/(2*kappa + 2))*G2;
  a1(isnan(a1)) = 0; % necessary because Inf*0 = NaN and we need 0
  
  b1 = a1.*Q1 + 1/2*G1 + ((Psi.^2)./(1 - Q1)).*G2;
  
  B1 = zeros(n1, K-1);
  for jj=1:K
      B1 = B1 + ((b1(:, jj) - (a1(:, jj).*Q1(:, jj)))*ones(1, K-1)).*UU1(:,:,jj);
  end
  
  if debug
      % constant terms in quadratic majorization
      c1 = a1.*(Q1.^2) + (1-kappa)/2*G1 + ...
          ((Psi.^2).*(1 + 2.*Q1./(1 - Q1))).*G2;
  
      TT1 = a1.*(Q1.^2) - 2*b1.*Q1 + c1;
      TT2 = H1;
      D = sum(sum(abs(TT1 - TT2)))/n1;
      if D > eps
          fprintf('\tSecond Case 1: diff = %15.16f\n', D);
      end
  
      % Case 1 constant terms in total majorization
      Gamma1 = 1/n * sum(P1.*sum(c1.*R1,2));
      Gamma1 = Gamma1 + 1/n * sum(P1.*sum(ZV1.^2,2).*sum(a1.*R1,2));
      Gamma1 = Gamma1 - 1/n * sum(P1.*sum(a1.*(Q1.^2).*R1,2));
  
      clear c1 TT1 TT2 D
  end
  
  %% Now create all matrices for the Case 2 objects.
  % We can now safely delete a number of matrices from memory
  clear G1 G2 G3 Phi Psi iPhi b1
  
  G1a = (Q2 <= (p+kappa-1)/(p-2));
  G2a = (Q2 <= 1)&(~G1a);
  G3a = ~(G1a|G2a);
  
  G1b = (Q2 <= -kappa);
  G2b = (Q2 <= 1) & (~G1b);
  G3b = ~(G1b|G2b);
  
  Phi = 1 - Q2 - (kappa+1)/2;
  Psi = (1 - Q2)/sqrt(2*kappa + 2);
  if p~=2
      Chi = (p*Q2 + kappa - 1)/(p - 2);
  end
  
  omega = (1/p)*(sum((H2.^p).*R2,2)).^(1/p - 1);
  
  if debug
      % First majorization test (p-th root)
      TT1 = sum((H2.^p).*R2,2).^(1/p);
      TT2 = omega.*sum((H2.^p).*R2,2) + (1 - 1/p)*(sum((H2.^p).*R2,2)).^(1/p);
      D = sum(abs(TT1 - TT2))/n2;
      if D > eps
          fprintf('\tFirst Case 2: diff = %15.16f\n', D);
      end
  end
  
  % Some parameters are different when p = 2, we recognize this here.
  if p~=2
      a2 = (1/4 * p^2 * Phi.^(p-2)).*G1a + ...
          (1/4 * p * (2*p - 1) * ((kappa+1)/2)^(p-2)).*G2a + ...
          (1/4 * p^2 * (p*Phi/(p-2)).^(p-2)).*G3a;
      a2(isnan(a2)) = 0; % We need Inf*0 = 0.
  else
      a2 = 3/2*ones(n2, K);
  end
  
  b2 = (a2.*Q2 + 1/2*p*(Phi.^(p-1))).*G1b + ...
      (a2.*Q2 + p*(Psi.^(2*p))./(1 - Q2)).*G2b;
  if p~=2
      b2 = b2 + (a2.*Chi + 1/2*p*(p*Phi/(p-2)).^(p-1)).*G3b;
  else
      b2 = b2 + 3/2*Q2.*G3b;
  end
  
  B2 = zeros(n2, K-1);
  for jj=1:K
      B2 = B2 + ((b2(:, jj) - (a2(:, jj).*Q2(:, jj)))*ones(1, K-1)).*UU2(:,:,jj);
  end
  
  if debug
      c2 = (a2.*(Q2.^2) + Phi.^p + p*Q2.*(Phi.^(p-1))).*G1b + ...
          (a2.*(Q2.^2) + (Psi.^(2*p)).*(1 + (2*p*Q2)./(1 - Q2))).*G2b;
      if p~=2
          c2 = c2 + (a2.*(Chi.^2) + p*Chi.*(p*Phi/(p - 2)).^(p-1) + (p*Phi/(p-2)).^p).*G3b;
      else
          c2 = c2 + 3/2*(Q2.^2).*G3b;
      end
  
      TT1 = a2.*(Q2.^2) - 2*b2.*Q2 + c2;
      TT2 = H2.^p;
      D = sum(sum(abs(TT1 - TT2)))/n2;
      if D>eps
          fprintf('\tSecond Case 2: diff = %15.16f\n', D);
      end
  
      % Case 2 constant terms in majorization
      Gamma2 = 1/n * (1 - 1/p) * sum(P2.*(sum((H2.^p).*R2,2).^(1/p)));
      Gamma2 = Gamma2 + 1/n * sum(P2.*omega.*sum(c2.*R2,2));
      Gamma2 = Gamma2 + 1/n * sum(P2.*omega.*sum(ZV2.^2,2).*sum(a2.*R2,2));
      Gamma2 = Gamma2 - 1/n * sum(P2.*omega.*sum((Q2.^2).*a2,2));
  end
  
  %% Collect the two classes in a single matrix and calculate update
  A(C1i, :) = P1.*sum(a1.*R1,2);
  A(C2i, :) = P2.*omega.*sum(a2.*R2,2);
  
  B(C1i, :) = (P1*ones(1, K-1)).*B1;
  B(C2i, :) = ((P2.*omega)*ones(1, K-1)).*B2;
  
  A = 1/n*A;
  B = 1/n*B;
  
  J = eye(m+1); J(1,1) = 0;
  
  ZAZ = Z'*((A*ones(1, m+1)).*Z);
  newV = (ZAZ + lambda*J)\(ZAZ*V + Z'*B);
  %% temp
  lhs = ZAZ + lambda*J;
  rhs = ZAZ*V + Z'*B;
  fprintf("lhs:\n");
  disp(lhs);
  fprintf("rhs:\n");
  disp(rhs);
  fprintf("newV:\n");
  disp(newV);
  %% endtemp
  
  if debug
      oldL = trace((V - 2*V)'*Z'*diag(A)*Z*V) - 2*trace(B'*Z*V) + Gamma1 + ...
          Gamma2 + lambda*trace(V'*J*V);
  else
      oldL = 0;
  end

end

function U = SimplexGen(K)
  U = zeros(K, K-1);
  for ii=1:K
      for jj=1:K-1
          if ii<=jj
              U(ii,jj) = -1/sqrt(2*(jj^2 + jj));
          elseif ii==jj+1
              U(ii,jj) = jj/sqrt(2*(jj^2 + jj));
          end
      end
  end
end