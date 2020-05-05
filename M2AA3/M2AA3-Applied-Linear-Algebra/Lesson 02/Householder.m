function [Q,R,O] = Householder(A)
% Compute full QR factorisation using householder algorithm.

[m, n] = size(A); % Calculate matrix column and rows
% Here m >= n.

Q = eye(m)

for i = 1:n
    col = A(i:m,i);
    ei = zeros(length(col),1);
    ei(1) = 1;
    v = col + mysign(col(1))*norm(col).*ei  % Choose Reflection Vector 
    out = (v*v')/(v' * v)
    Pred = eye(length(v)) - 2*out           % Reflection Matrix
    P = blkdiag(eye(i-1), Pred)
    A = P*A
    Q = Q*P
end
R = A;
O = Q'*Q;
end

function y = mysign(x)
    if x >= 0   % As a convention, sign(0) = 0
        y = 1;
    else
        y = -1;
    end
end

