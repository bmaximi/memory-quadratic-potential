function gammaj = solve_constr_lsq(lambdaj,y,t,Omega,opt)
    
    % Set necessary parameters
    if ~isfield(opt,'acf')
        opt.acf = 'v';
    end

    dim = size(y,1);
    m = length(lambdaj);
    I = eye(dim*dim);
    orig = reshape(1:dim^2,dim,dim);
    permut = orig';
    P = I(permut(:),:);
    
    upperT = triu(orig,1);
    evenIndex = upperT(upperT~=0);
    
    b = y(:);
    
    % constr eq. matrices
    Y0 = kron(ones(1,m),I);
    Y1 = kron(lambdaj.',I);
    Y2 = 1i*kron(lambdaj.^2.',I-P);
    Y3 = kron(lambdaj.^3.',I+P);
    Psi0 = -kron(lambdaj.^(-1).',I);
    Psi1 = kron(lambdaj.^(-2).',I);
    Psi2 = kron(lambdaj.^(-3).',I+P);
    
    % switch 
    if opt.acf == 'p'
        G = -kron( (lambdaj(:).^(-2).*exp(lambdaj(:) .* t)) .', I );
    else
        G = kron( exp(lambdaj(:) .* t) .', I );
    end
    
    B = [Psi0; Y0; Y1; Y2(evenIndex,:);];
    d = [zeros(dim^2,1); reshape(eye(dim),[],1); zeros(dim^2+length(evenIndex),1)];

    if Omega~=inf & size(B,1) < dim*dim*m
        B = [Psi1; B];
        d = [-reshape(inv(Omega),[],1);d];
    end

    a = constr_lsq(G,b,B,d);

    % check if sdp constraints are satisfied
    sdp1 = min(eig(reshape(Y3*a,dim,dim))) > 0;
    sdp2 = min(eig(reshape(Psi2*a,dim,dim))) > 0;
    
    counter = 1; % counter to avoid dead loop below
    
    % if sdp constraints are violated, add it to list of constraints and
    % solve augmented optimization problem
    while (~sdp1 || ~sdp2) & counter <= 2
        cvx_begin quiet sdp
            variable a(m*dim*dim) complex;
            minimize( norm(b-G*a, 'fro') );
            subject to
            B*a == d;
            if ~sdp1
                reshape(Y3*a,dim,dim) == semidefinite(dim);
            end
            if ~sdp2
                reshape(Psi2*a,dim,dim) == semidefinite(dim);
            end
        cvx_end

        if convertCharsToStrings(cvx_status) == 'Infeasible'
            error('*** Semidefinite Problem is infeasible');
        end

        sdp1 = min(eig(reshape(Y3*a,dim,dim))) > 0;
        sdp2 = min(eig(reshape(Psi2*a,dim,dim))) > 0;
        counter = counter + 1;
    end
    
    if ~sdp1 || ~sdp2
        cvx_begin quiet sdp
            variable a(m*dim*dim) complex;
            minimize( norm(b-G*a) );
            subject to
            B*a == d;
            reshape(Y3*a,dim,dim) == semidefinite(dim);
            reshape(Psi2*a,dim,dim) == semidefinite(dim);
        cvx_end

        if convertCharsToStrings(cvx_status) == 'Infeasible'
            error('*** Semidefinite Problem is infeasible');
        end
    end
    gammaj = reshape(a,dim,dim,[]);
    
end

function x = constr_lsq(A,b,B,d)
% solves the constrained least-squares problem
%
% minimize ||Ax - b||_2     subject to Bx = d
%
% Input: A,b:           Parameters of Least-Squares-problem
%        B,d:           Parameters of constraints
%
% Output: x:            Least-Square solution

% Solve constrained LSQ-Problem with nullspace-method

% check for unique solvability of lsq problem
S = svd([A;B]);
if sum(S(S<1e-14)) ~= 0
    warning('Constrained LSQ has no unique minimizer');
end

p = size(B,1);
n = size(B,2);
[QB,RB] = qr(B');
Q1 = QB(:,1:p);
Q2 = QB(:,p+1:end);
x1 = Q1*(RB(1:p,1:p)'\d);
[QA,RA] = qr(A*Q2);
c = QA'*(b-A*x1);
y2 = RA(1:n-p,1:n-p)\c(1:n-p);

x = x1 + Q2*y2;

end

