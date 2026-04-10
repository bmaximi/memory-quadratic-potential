function X=Riccati(P,G,Q)
%
%   X = Riccati(A,G,Q) solves the algebraic Riccati equation of the form:
%       P*X + X*P' + X*G*X + Q = 0, where X is symmetric.
%
n=size(P,1);
Z=[P' G; -Q -P];
[U1,S1]=schur(Z);
[U,~]=ordschur(U1,S1,'lhp');
X=U(n+1:end,1:n)/U(1:n,1:n);