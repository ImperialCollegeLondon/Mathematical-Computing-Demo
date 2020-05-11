function [x,r] = Solve_LUP(A,b,haveP)
% Solve Ax=b via LU if haveP == 0 or LUP if haveP == 1.

[m,n] = size(A);
[L,U,P] = LUP(A,haveP)
bnew = P*b

% Forward substitution phase
x = zeros(n,1);          % Pre-allocation for x
x(1) = bnew(1)/L(1,1)
for i = 2:n
    v = dot(x,L(i,:))
    x(i) = (bnew(i)-v)/L(i,i)
end

% Back substitution phase
bnew = x
x = zeros(n,1);
x(n) = bnew(n)/U(n,n)
for i = 1:n-1
    v = dot(x,U(n-i,:))
    x(n-i) = (bnew(n-i)-v)/U(n-i,n-i)
end

% Obtaining residue
rvec = A*x - b
r = norm(rvec)
end

