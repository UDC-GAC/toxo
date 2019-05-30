function A = nfold(s)
%NFOLD N-times cartesian product of a vector with itself.
%   Source code provided by John D'Errico at:
%   https://mathworks.com/matlabcentral/answers/282777-n-fold-cartesian-product#answer_220945

n = length(s);
[A{n:-1:1}] = ndgrid(s{n:-1:1});
A = reshape(cat(n, A{:}), [], n);
end
