%%
% adj_to_trans.m
%
% Given a (square) adjacency matrix A, converts A to a transition matrix 
% using the softmax function on connected edges.
% Returns transition matrix T.

function [T] = adj_to_trans(A)

nrow = length(A); % Number of rows
T = A;

for i = 1:nrow
    row = A(i,:);
    k = find(row); % Find the index which value is not 0.
    T(i,k) = softmax(row(k).').';
end
end

