function [p] = genotype_probabilities(maf, order)
% GENOTYPE_PROBABILITIES Compute the probabilities associated with each genotype for a given order and MAF.
m = sym(maf);
M = 1 - m;
p = prod(toxo.nfold([M^2, 2 * M * m, m^2], order), 2);
end

