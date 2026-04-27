function out = enorm(X);

% this function takes the elements of the column vector X and produces the
% norm of each element in the column vector out.

nx = zeros(size(X,1),1)*NaN;

for i = 1:size(X,1)
    nx(i) = norm(X(i,1));
end
out = nx;