function A = nfold(vector, n)
%NFOLD Summary of this function goes here
%   https://www.mathworks.com/matlabcentral/answers/282777-n-fold-cartesian-product
s = repmat({vector}, 1, n);
[A{n:-1:1}] = ndgrid(s{:});
A = reshape(cat(n, A{:}), [], n);
end
