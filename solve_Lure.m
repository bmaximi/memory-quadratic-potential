function [Sigma, G] = solve_Lure(A,d)

% Partition of A
B = A(1:d,d+1:end)';
C = -A(d+1:end,1:d);
A0 = A(d+1:end,d+1:end);

N = size(A0,1);

% Reduction step
X = basisTrafo(B,C);
XAX = X*A0/X;
D1 = -XAX(1:d,1:d);
B1 = XAX(1:d,d+1:end)';
C1 = -XAX(d+1:end,1:d);
A1 = XAX(d+1:end,d+1:end);

R1 = D1 + D1';

if min(eig(R1))<0
end

if size(A1,1) > 0 && max(real(eig(A1)))>=0
end

% solve Riccati equation
P = A1 - C1*(R1\B1');
Q1 = B1*(R1\B1');
Q2 = C1*(R1\C1');

Sigma1 = Riccati(P,Q1,Q2);

err = norm(P*Sigma1+Sigma1*P' + Sigma1*Q1*Sigma1 + Q2)/norm(P*Sigma1);

if err > 1e-6
end
if size(Sigma1,1) > 0 && norm(Sigma1-Sigma1') > 1e-5
end

% transformation to full system
K1 = chol(R1)';
L1 = (C1 - Sigma1*B1)/K1';

XL = [K1; L1]; L  = X\XL;
XSigma0X = [eye(d), zeros(d,N-d); zeros(N-d,d), Sigma1];
Sigma0 = X\XSigma0X/X';
Sigma = (Sigma0+Sigma0')/2;

G = [zeros(d); L];

end