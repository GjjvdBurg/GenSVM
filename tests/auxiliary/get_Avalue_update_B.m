clear;
more off;
rand('state', 21831);

n = 5;
K = 3;
m = 3;

y = [2 1 3 2 3]';

R = zeros(n, K);
I = eye(K);
for i=1:n
  R(i, :) = I(y(i, :), :);
end
R = ~logical(R);

omega = [
  0.7394076262220608
  0.7294526264247443
  0.6802499471888741
  0.6886792032441273
  0.8695329737474283
];


H = zeros(5, 3);

H(1, 1) = 0.8465725800087526;
H(1, 2) = 1.2876921677680249;
H(1, 3) = 1.0338561593991831;
H(2, 1) = 1.1891038526621391;
H(2, 2) = 0.4034192031226095;
H(2, 3) = 0.0;
H(3, 1) = 0.5;
H(3, 2) = 0.0;
H(3, 3) = 1.1;
H(4, 1) = 0.0;
H(4, 2) = 0.0;
H(4, 3) = 0.8189336598535132;
H(5, 1) = 0.6164203008844277;
H(5, 2) = 0.2456444285093894;
H(5, 3) = 0.8184193969741095;