%%
% min_k.m
%

function [vals, idxs] = min_k(vec, k)
    k = min(k, length(vec));
    [sorted_vals, sorted_idxs] = sort(vec);
    vals = sorted_vals(1:k);
    idxs = sorted_idxs(1:k);
end

