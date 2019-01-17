function [p] = genotype_probabilities(maf, order)
% GENOTYPE_PROBABILITIES Compute the probabilities associated with each genotype for a given order and MAF.

p = prod(toxo.nfold([(1 - maf)^2, 2 * maf * (1 - maf), maf^2], order), 2);
end

