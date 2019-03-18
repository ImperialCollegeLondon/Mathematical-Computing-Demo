function [Q,R,O] = mgs(A)
% Input: Matrix A
% Output:
% Matrix Q - The orthogonal matrix which contains the transformed vectors
% Matrix R - An upper-right triangle matrix
% Matrix O - for Checking, which is Q' * Q - I(n)

% Initialization
[m, n] = size(A); % Calculate matrix column and rows
V = A;            % Temporary Matrix
Q = zeros(m,n);   % Pre-allocation for Q

% MGS
% Loop for V (for first n column)
for i = 1:n-1
    Q(:,i) = V(:,i)/norm(V(:,i));
    for j = i+1 : n
        V(:,j) = V(:,j) - dot(V(:,j),Q(:,i))*Q(:,i);
    end
end

% Normalizing Vectors
for i = 1:n
    Q(:,i) = V(:,i)/norm(V(:,i));
end

% Loop for R (for first n column)
R = zeros(n);
for q = 1:n %column
    for p = 1:q %row
        R(p,q) = dot(Q(:,p),A(:,q));
    end
end

% Checking
O = (Q')*Q - eye(n);
end