function b = ctx_dec2col(X, n)
if nargin < 2
    n = max(X(:));
end
b = zeros(length(X),n); % Create an empty array w/ as many rows as X, and as many cols as the max val of X

for i=1:length(X) 
    b(i,X(i)) = 1;   % In that array, for each row and the corresponding column where X=that row, put a 1
end
end
