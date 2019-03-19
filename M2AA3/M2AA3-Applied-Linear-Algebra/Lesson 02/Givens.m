function [Q,R,O] = Givens(A)
% Input: Matrix A
% Output:
% Matrix Q - The orthogonal matrix which contains the transformed vectors
% Matrix R - An upper-right triangle matrix
% Matrix O - for Checking, which is Q' * Q - I(n)

% Initialization
[m, n] = size(A); % Calculate matrix column and rows
G = eye(m);       % Initialization of G
R = A;            % Initializing R

% Loop - Column then Row
for p = 1:n                 % Column of entries to reduce to 0
    for q = p+1:m           % Row of entries to reduce to 0
        Gpq = eye(m); % Constructing Gpq
        x = R(p,p);
        y = R(q,p);
        c = x/sqrt(x^2 + y^2);
        s = -y/sqrt(x^2 + y^2);
        Gpq(p,p) = c;
        Gpq(q,q) = c;
        Gpq(p,q) = -s;
        Gpq(q,p) = s;
        G = Gpq * G;
        R = Gpq * R;
    end
end
Q = G';
O = Q' * Q - eye(m);
end