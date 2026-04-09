function poles = calcPoles(W,z)
    
    d = size(W,2);
    m = length(z);
    
    B = eye(m*d+d); B(1:d,1:d) = 0; E = 0*B;
    if size(W,1)/d == m
        E(1:d,d+1:end) = W.';
    else
        E(1:d,:) = W.';
    end
    E(d+1:end,1:d) = kron(ones(m,1),speye(d));
    E(d+1:end,d+1:end) = kron(diag(z),speye(d));
    
    poles = eig(E,B);
    poles = poles(abs(poles)<1);
end
        
    