%%
% get_lap.m
%
% Compute the graph Laplacians
% Input: adjacency matrix A


function [L] = get_lap(A)
    d = size(A, 1);

    % Compute node degrees.
    degrees = sum(A, 2);
    
    % Compute graph Laplacians.
    D = diag(degrees);
    L = D - A;
end

   