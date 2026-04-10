function out = markovianAppPot(y,t,g,r,Omega,opt)
    
% this function computes a Markovian approximation to a given grid using
% rational approximation and (semidefinite) optimization
%
% Input:    y,t: 2n data points and time grid with width tau
%           g,r: number of grid points and radius for generating function
%           Omega: stiffness matrix of harmonic potential (inf if not known)
%           opt: contains optional parameters like aaa tolerance and used
%           acf data (standards are aaatol = 1e-4 and acf = 'v')

% Output:   out: struct containing the Prony series and Markov parameters

tau = t(2);
dim = size(y,1);

% generating function + rational approximation for poles
[omega, F] = genFunc(y,r,g);
out.poles = blockaaa(omega,F,opt);

% exponents
out.lambdaj = calcExponents(out.poles,tau);

% coefficients with optimization problem
out.gammaj = solve_constr_lsq(out.lambdaj,y,t,Omega,opt);

% set up system matrix
out.A = constrSysmat(out.lambdaj,out.gammaj);
m = size(out.A,1);

% solve Lur'e system
[Sigma0, G0, out.flag] = solve_Lure(out.A(1:m-dim,1:m-dim),dim);

out.Sigma = [eye(dim), zeros(dim,m-dim);...
    zeros(m-2*dim,dim), Sigma0, zeros(m-2*dim,dim);...
    zeros(dim,m-dim), -inv(out.A(1:dim,m-dim+1:m))];
out.G = [G0; zeros(dim)];

if out.flag == 0
    fprintf('Reconstruction successful\n')
end

end