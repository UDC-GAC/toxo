function [p] = genotype_probabilities(mafs)
% GENOTYPE_PROBABILITIES Compute the probabilities associated with all genotype combinations given each MAF.

m = sym(mafs);
M = 1 - m;
s = arrayfun(@(m, M) {[M^2, 2 * M * m, m^2]}, m, M);
p = prod(toxo.nfold(s), 2);
end
