function poles = blockaaa(Z,F,opt)
% Computes a rational approximation to a given set of data using the AAA
% method.

    % Set necessary parameters
    if ~isfield(opt,'aaatol')
        opt.aaatol = 1e-4;
    end

    % init aux variables
    Z = Z(:);
    M = length(Z);
    J = 1:M;
    d = size(F,1);

    % initial iteration with real data points
    z = [Z(1); Z(M/2+1)];
    f(:,:,1) = F(:,:,1); f(:,:,2) = F(:,:,M/2+1);
    J(J==1) = []; J(J==M/2+1) = [];
    m=2;

    % optimize weights
    L = constrLoewner(J,F,f,Z,z);
    [~,~,V] = svd(L,0);

    w = symmWeights(V);

    % compute rat app
    R = computeRatApp(F,f,J,w,Z,z);

    errvec = 0*Z;
    for i=1:M
        errvec(i) = norm(F(:,:,i) - R(:,:,i),'fro');
    end
    [~,ind] = max(errvec);
    err = inf;

    while (m < M/2 && err > opt.aaatol)
        % choose next support point
        m = m+2;
        cc = M-ind+2;
        z = [z; Z(ind); Z(cc)];
        f(:,:,m-1) = F(:,:,ind); f(:,:,m) = F(:,:,cc);
        J(J==ind) = []; J(J==cc) = [];

        L = constrLoewner(J,F,f,Z,z);
        [~,~,V] = svd(L,0);

        w = symmWeights(V);
        R = computeRatApp(F,f,J,w,Z,z);

        for i=1:M
            errvec(i) = norm(F(:,:,i) - R(:,:,i),'fro');
        end
        [err, ind] = max(errvec);
    end

    poles = calcPoles(w,z);
end

function L = constrLoewner(J,F,f,Z,z)
    d = size(F,1);
    m = length(z);

    L = zeros(length(J)*d^2,m);
    for j=1:m
        tmp = zeros(length(J)*d,d);
        for k=1:length(J)
            tmp((k-1)*d+1:k*d,:) = ( F(:,:,J(k)) - f(:,:,j) ) / (Z(J(k)) - z(j));
        end
        L(:,j) = tmp(:);
    end 
end

function symmW = symmWeights(V)
    n = size(V,1);
    w1 = V(:,end);
    w2 = w1;
    for i=3:2:n
        w2(i) = w1(i+1);
        w2(i+1) = w1(i);
    end
    symmW = w1 + conj(w2);
    symmW = symmW / norm(symmW,'fro');
end

function R = computeRatApp(F,f,J,w,Z,z)
    R = F;
    d = size(F,1);

    for k=1:length(J)
        N = zeros(d); D = 0;
        for l=1:length(z)
            N = N +(w(l) * f(:,:,l)) / (Z(J(k)) - z(l));
            D = D + w(l) / (Z(J(k)) - z(l));
        end
        R(:,:,J(k)) = N / D;
    end
end

