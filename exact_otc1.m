%%
% exact_otc1.m
%
% Run ExactOTC.

function [exp_cost, P, stat_dist] = exact_otc1(Px, Py, c)
dx = size(Px, 1);
dy = size(Py, 1);

P_old = ones(dx*dy);
P = get_ind_tc(Px, Py);
iter_ctr = 0;
while max(max(abs(P-P_old))) > 1e-10
    iter_ctr = iter_ctr + 1;
    P_old = P;
    
    % Transition coupling evaluation.
    [g, h] = exact_tce(P, c);
    
    % Transition coupling improvement.
    P = exact_tci(g, h, P_old, Px, Py);
    
    % Check for convergence.
    if all(all(P == P_old))
        stat_dist = get_stat_dist(P);
        stat_dist = reshape(stat_dist, dy, dx)';
        exp_cost = sum(stat_dist .* c, 'all');
        return
    end
end
end
