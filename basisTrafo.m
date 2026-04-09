function X = basisTrafo(U,V)
    % Constructs a transformation matrix such that
    % XV = ES^(1/2) and X^-T E = US^(-1/2)
    %
    % Input:    U,V: matrix pair with S = U'*V being symmetric and pos.
    %                def.
    % Output:   X: transformation matrix that satisfies above identity

    S = U'*V;
    [Q,~,~] = svd(V);
    d = size(S,1);
    
    X = [S^(-1/2)*U'; Q(:,d+1:end)'];
end