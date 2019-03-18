function [Q,R,O] = cgs(A)
% Input: Matrix A
% Output:
% Matrix Q - The orthogonal matrix which contains the transformed vectors
% Matrix R - An upper-right triangle matrix
% Matrix O - for Checking, which is Q' * Q - I(n)

% Initialization
[m, n] = size(A); % Calculate matrix column and rows
V = zeros(m,n);   % Preallocation for V, matrix for non-normalized vectors.
Q = zeros(m,n);   % Pre-allocation for Q
V(:,1) = A(:,1);
Q(:,1) = A(:,1)/norm(A(:,1));  %Obtaining q1

% CGS
for j = 2:n
    
    % Summation Loop
    vec = A(:,j);
    for i = 1:j-1
        vec = vec - dot(A(:,j),Q(:,i))*Q(:,i);
    end
    
    % Storing vectors in V and Q
    V(:,j) = vec;
    Q(:,j) = vec / norm(vec);
end
   
% Loop for R
R = zeros(n);     % Pre-allocation for Q
for q = 1:n       % column index
    for p = 1:q-1 % row index
        R(p,q) = dot(Q(:,p),A(:,q));  % Non Diagonal Entries
    end
    R(q,q) = norm(V(:,q));            % Diagonal Entries
end

% Checking
O = (Q')*Q - eye(n);
end