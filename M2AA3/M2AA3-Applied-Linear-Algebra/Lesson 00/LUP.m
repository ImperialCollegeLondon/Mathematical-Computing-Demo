function [L,U,P] = LUP(A,haveP)
% Perform LUP Decomposition if haveP is 'true', otherwise perform LU
% Decomposition. Only valid for square matrix A.

[m,n] = size(A);
L = eye(m);
U = A;
P = eye(m);

if haveP == 1
    for k = 1:m-1
        % pivoting
        [val, swapidx] = max(abs(U(k:end,k)));
        swapidx = swapidx + k - 1
        U([k swapidx],k:end) = U([swapidx k],k:end)
        L([k swapidx],1:k-1) = L([swapidx k],1:k-1)
        P([k swapidx],:) = P([swapidx k],:)
        for j = k+1:m
            L(j,k) = U(j,k)/U(k,k)
            U(j,:) = U(j,:) - L(j,k)*U(k,:)
        end
    end
else
    for k = 1:m-1
        for j = k+1:m
            L(j,k) = U(j,k)/U(k,k)
            U(j,:) = U(j,:) - L(j,k)*U(k,:)
        end
    end
end