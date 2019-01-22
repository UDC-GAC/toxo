function A = nfold(vector, n)
%NFOLD N-times cartesian product of a vector with itself.
%   Source code provided by John D'Errico at:
%   https://mathworks.com/matlabcentral/answers/282777-n-fold-cartesian-product#answer_220945

s = repmat({vector}, 1, n);
[A{n:-1:1}] = ndgrid(s{:});
A = reshape(cat(n, A{:}), [], n);
end
