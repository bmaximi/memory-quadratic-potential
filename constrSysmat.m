function smat = constrSysmat(lambdaj,gammaj)

    d = size(gammaj,1);
    m = length(lambdaj);
    r = 1;
    J = []; USigma = []; VSigma = [];

    % real exponents
    while abs(imag(lambdaj(r))) < 1e-8
        [U,S,V] = svd(gammaj(:,:,r)); S = diag(S);
        for i=1:d
            if S(i) < 1e-12
                break
            end
            J(end+1,end+1) = real(lambdaj(r));
            USigma(:,end+1) = sqrt(S(i))*real(U(:,i));
            VSigma(end+1,:) = sqrt(S(i))*real(V(:,i))';
        end
        r = r+1;
        if r > m
            break
        end
    end

    % complex exponents
    while r<m
        [U,S,V] = svd(gammaj(:,:,r)); S = diag(S);
        for i=1:d
            if S(i) < 1e-12
                break
            end
            J(end+1:end+2,end+1:end+2) = [real(lambdaj(r)), imag(lambdaj(r)); -imag(lambdaj(r)), real(lambdaj(r))];
            USigma(:,end+1) = sqrt(2*S(i))*real(U(:,i));
            USigma(:,end+1) = sqrt(2*S(i))*imag(U(:,i));
            VSigma(end+1,:) = sqrt(2*S(i))*real(V(:,i)');
            VSigma(end+1,:) = -sqrt(2*S(i))*imag(V(:,i)');
        end
        r = r+2;
    end

    X = basisTrafo(USigma',VSigma);
    A = X*J*inv(X);

    % bring bottom right block into Jordan canonical form
    [U,D] = eig(A(d+1:end,d+1:end),'vector');
    [~,ind] = sort(abs(real(D)),'descend');
    U = U(:,ind); D = diag(D(ind));
    
    % set up trafo matrix
    T = 0.5*[1 -1i; 1 1i];
    X = zeros(size(D,1)); k=1;
    while k <= size(D,1)
        if abs(imag(D(k,k))) < 1e-10
            X(k,k) = 1;
            k = k+1;
        else
            X(k:k+1,k:k+1) = T;
            k = k+2;
        end
    end
    
    % apply transformation
    smat = [zeros(d), A(1:d,d+1:end)*U*X; ...
    (U*X)\A(d+1:end,1:d), (U*X)\A(d+1:end,d+1:end)*U*X];
    
    smat(1:d,end-d+1:end) = smat(1:d,end-d+1:end)*smat(end-d+1:end,1:d);
    smat(end-d+1:end,1:d) = eye(d);

end