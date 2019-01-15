
function test_get_alpha_beta()

  n = 8;
  K = 3;
  m = 3;

  kappa = 0.5;
  p = 1.1;

  rho = ones(n, 1);

  V = zeros(m+1, K-1);

  V(1, 1) = -0.7593642121025029;
  V(1, 2) = -0.5497320698504756;
  V(2, 1) = 0.2982680646268177;
  V(2, 2) = -0.2491408622891925;
  V(3, 1) = -0.3118572761092807;
  V(3, 2) = 0.5461219445756100;
  V(4, 1) = -0.3198994238626641;
  V(4, 2) = 0.7134997072555367;

  Z = zeros(n, m+1);

  Z(1, 1) = 1.0000000000000000;
  Z(1, 2) = 0.6437306339619082;
  Z(1, 3) = -0.3276778319121999;
  Z(1, 4) = 0.1564053473463392;
  Z(2, 1) = 1.0000000000000000;
  Z(2, 2) = -0.8683091763200105;
  Z(2, 3) = -0.6910830836015162;
  Z(2, 4) = -0.9675430665130734;
  Z(3, 1) = 1.0000000000000000;
  Z(3, 2) = -0.5024888699077029;
  Z(3, 3) = -0.9649738292750712;
  Z(3, 4) = 0.0776560791351473;
  Z(4, 1) = 1.0000000000000000;
  Z(4, 2) = 0.8206429991392579;
  Z(4, 3) = -0.7255681388968501;
  Z(4, 4) = -0.9475952272877165;
  Z(5, 1) = 1.0000000000000000;
  Z(5, 2) = 0.3426050950418613;
  Z(5, 3) = -0.5340602451864306;
  Z(5, 4) = -0.7159704241662815;
  Z(6, 1) = 1.0000000000000000;
  Z(6, 2) = -0.3077314049206620;
  Z(6, 3) = 0.1141288036288195;
  Z(6, 4) = -0.7060114827535847;
  Z(7, 1) = 1.0000000000000000;
  Z(7, 2) = 0.6301294373610109;
  Z(7, 3) = -0.9983027363627769;
  Z(7, 4) = -0.9365684178444004;
  Z(8, 1) = 1.0000000000000000;
  Z(8, 2) = -0.0665379368401439;
  Z(8, 3) = -0.1781385556871763;
  Z(8, 4) = -0.7292593770500276;

  y = [ 2 , 1 , 3 , 2 , 3 , 3 , 1 , 2 ];

  % Calculate ZV
  ZV = Z*V;

  U = SimplexGen(3);

  % Calculate errors

  Q = zeros(n, K);
  for i=1:n
    for j=1:K
      if (j == y(i))
        continue
      end

      Q(i, j) = ZV(i, :) * (U(y(i), :) - U(j, :))';
    end
  end

  % Calculate huber

  H = zeros(n, K);
  for i=1:n
    for j=1:K
      q = Q(i, j);
      value = 0.0;
      if (q <= -kappa)
        value = 1 - q - (kappa + 1) / 2;
      elseif (q <= 1.0)
        value = 1/(2*kappa + 2) * (1 - q)*(1 - q);
      end
      H(i, j) = value;
    end
  end

  % The index we're computing alpha / beta for
  % Note the offset for 0-based and 1-based
  i = 1;

  is_simple = true;
  value = 0;
  for j=1:K
    if (j == y(i))
      continue
    end
    h = H(i, j);
    value += (h > 0);
    if (value > 1)
      is_simple = false;
      break;
    end
  end
  fprintf('Is simple: %d\n', is_simple);

  omega = 0.0;
  for j=1:K
    if (j == y(i))
      continue
    end
    h = H(i, j);
    omega = omega + h^p;
  end
  omega = (1.0/p) * (omega ^ (1.0/p - 1.0));

  if (is_simple)
    omega = 1.0;
  end

  fprintf('Omega: %.16f\n', omega);

  thebeta = zeros(1, K-1);
  thealpha = 0.0;

  for j=1:K
    if (j == y(i))
      continue
    end

    if (is_simple)
      q = Q(i, j);
      if (q <= -kappa)
        a = 0.25 / (0.5 - kappa / 2.0 - q);
        b_aq = 0.5;
      elseif (q <= 1.0)
        a = 1.0/(2.0*kappa + 2.0);
        b_aq = (1 - q)*(a);
      else
        a = -0.25 / (0.5 - kappa/2.0 - q);
        b_aq = 0;
      end
    else
      q = Q(i, j);
      a2g2 = 0.25 * p * (2.0 * p - 1.0) * ((kappa + 1)/2)^(p - 2);

      if (2 - p < 1e-2)
        if (q <= -kappa)
          b_aq = 0.5 - kappa/2.0 - q;
        elseif (q <= 1)
          b_aq = ((1.0 - q)^3.0) / (2.0 * (kappa + 1.0)^2.0);
        else
          b_aq = 0;
        end
        a = 1.5;
      else
        if (q <= (p + kappa - 1)/(p - 2))
          a = 0.25 * p^2 * (0.5 - kappa/2.0 - q)^(p - 2)
        elseif (q <= 1.0)
          a = a2g2;
        else
          a = 0.25 * p^2.0 * ((p/(p - 2))*(0.5 - kappa/2.0 - q))^(p - 2);
          b_aq = a*(2.0*q + kappa - 1.0)/(p - 2) + 0.5*p*((p/(p - 2)*(0.5 - kappa/2 - q)))^(p - 1);
        end
        if (q <= -kappa)
          b_aq = 0.5 * p * (0.5 - kappa/2.0 - q)^(p - 1);
        elseif (q <= 1.0)
          b_aq = p * ((1 - q)^(2*p - 1.0))/((2*kappa + 2)^p);
        end
      end
    end

    b_aq = b_aq * rho(i) * omega * (1 / n);

    uu_row = U(y(i), :) - U(j, :);

    thebeta = b_aq * uu_row + thebeta

    thealpha = thealpha + a;

  end

  thealpha = thealpha * (omega * rho(i) * 1.0 / n);

  fprintf('alpha:\n%.16f\n\n', thealpha);
  fprintf('beta:\n');
  disp(thebeta);

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
