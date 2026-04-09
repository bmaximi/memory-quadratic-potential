function [omega, F] = genFunc(y,r,g)
    % Evaluates the generating function on an equiangular grid 
    % f(z) = sum_n=0^N y_n z^(-n-1)
    %
    % Input: r,g        radius and number of grid points
    %        y          data points
    %
    % Output: omega, F  Grid and data points of generating function
    
    d = size(y,1);
    N = size(y,3);

    omega = r*exp(2*pi*1i*(0:g-1)/g);
    F = zeros(d,d,g);
    for i=1:g
        for j=1:N
            F(:,:,i) = F(:,:,i) + y(:,:,j)*omega(i)^(-j);
        end
    end
end