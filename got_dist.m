%%
% got_dist.m
%
% GOT distance between graphs with same number of nodes
% Input: Two Laplacian matrix
% Output: Wasserstein distance

function [W] = got_dist(A, B)
    
    % Check if A and B have same number of nodes
    if size(A, 1) ~= size(B, 1)
        error('Two graphs must have same number of nodes.')
    end 
    
    n = size(A, 1);
    
    % Adding 1 to zero eigenvalue
    % Does not change results, but is faster and more stable
    A = A + ones(n)/n;
    B = B + ones(n)/n;
    
    A_inv = pinv(A);
    B_inv = pinv(B);
    
    Root_A = sqrtm(A_inv);
 
    W = trace(A_inv) + trace(B_inv) - 2*trace(sqrtm(Root_A * B_inv * Root_A));
    
end

