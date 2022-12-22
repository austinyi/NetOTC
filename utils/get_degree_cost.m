%%
% get_degree_cost.m
%
% Compute the cost matrix of the squared difference in node degrees between
% graphs with weighted adjacency matrices D1 and D2.

function cost_mat = get_degree_cost(D1, D2)
    d1 = size(D1, 1);
    d2 = size(D2, 1);
    % Compute node degrees.
    degrees1 = sum(D1,2);
    degrees2 = sum(D2,2);
    % Construct matrix.
    cost_mat = zeros(d1, d2);
    for i=1:d1
        for j=1:d2
            cost_mat(i, j) = (degrees1(i) - degrees2(j))^2;
        end
    end
end