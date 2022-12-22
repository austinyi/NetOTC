%%
% get_best_stat_dist.m
%
% Compute best stationary distribution of a transition matrix given
% a cost matrix c. 
%
% Inputs:
% -P: transition matrix in R^(n x n)
% -c: cost matrix in R^(n x n)
%
% Outputs:
% -stat_dist: vector in R^n corresponding to best stationary distribution
% of P with respect to c.
% -exp_cost: the expected cost of stat_dist with respect to c, computed by
% simply taking the inner product of stat_dist and c.

function [stat_dist, exp_cost] = get_best_stat_dist(P, c)
    % Set up constraints.
    n = size(P,1);
    c = reshape(c', n, []);
    Aeq = [P' - eye(n);
           ones(1, n)];
    beq = [zeros(n, 1);
           1];
    lb = zeros(n,1);
    % Solve linear program.
    options = optimset('Display','off', 'TolCon', 1e-3, 'TolFun', 1e-6);
    options.Preprocess = 'none';
    %options = optimset('Display','off');
    [stat_dist, exp_cost] = linprog(c, [], [], Aeq, beq, lb, [], options);
    
    % In case the solver fails due to numerical underflow, try with
    % rescaling.
    alpha = 1;
    while isempty(stat_dist) & alpha >= 1e-10
        alpha = alpha/10;
        [stat_dist, exp_cost] = linprog(c, [], [], alpha*Aeq, alpha*beq, lb, [], options);
    end
    if isempty(stat_dist)
        error('Failed to compute stationary distribution.');
    end
end