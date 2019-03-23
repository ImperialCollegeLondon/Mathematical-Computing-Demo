function [L,O] = Cholesky(A)
% Input - A spd matrix
% Output: L and O = A - LL'

% Initialization
[m, n] = size(A); % Calculate matrix column and rows
L = [];     % Initializing L
R = A;

% Loop
for i = 1:n
    l = R(:,i)/sqrt(R(i,i));
    L = [L l];
    R = R - l*l';
end

O = A - L*L';
end

