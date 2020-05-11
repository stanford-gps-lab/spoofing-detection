function binPerm = binaryPermutations(n, k)
%ssPerm = subsetPermutations(n, k)
%   Generate matrix of all possible combinations choosing k out of n
%   options. Each row corresponds to an nchoosek permutation. Chosen values
%   are equal to 1, others equal to 0.

ii = nchoosek(1:n, k);
m = size(ii, 1);
binPerm = false(m, n);
binPerm(sub2ind([m, n], (1:m)'*ones(1, k), ii)) = 1;

end