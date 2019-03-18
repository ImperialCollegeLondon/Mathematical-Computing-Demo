function [x, r] = Solve_cgs(A, b)
% Input: Matrix A
% Output:
% Vector x - Solution Vector
% Matrix r - Residual Vector

% Initialization
[m, n] = size(A);
[Q, R, O] = cgs(A);      % Obtain QR factorization for A
bnew = Q' * b;
x = zeros(n,1);          % Pre-allocation for x
x(n) = bnew(n) / R(n,n); % Start from last row

for i = 1:n-1   % Row Index
    
    % Summation Loop
    v = 0;
    for j = 1:i % Column Index
        v = v + R(n-i,n-i+j)*x(n-i+j);
    end
    
    %Obtaining Entry
    x(n-i) = (bnew(n-i) - v)/R(n-i,n-i);

end

% Obtaining Residue
r = A*x - b;
rnorm = norm(r)
end

