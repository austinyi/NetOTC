%%
% entropic_otc.m
%
% Run EntropicOTC.

function [exp_cost, P, stat_dist] = entropic_otc(Px, Py, c, L, T, xi, sink_iter, get_sd)
dx = size(Px, 1);
dy = size(Py, 1);
max_c = max(max(c));
tol = 1e-5*max_c;

g_old = max_c*ones(dx*dy, 1);
g = g_old-10*tol;
P = get_ind_tc(Px, Py);
iter_ctr = 0;
while g_old(1) - g(1) > tol
    iter_ctr = iter_ctr + 1;
    P_old = P;
    g_old = g;
    
    % Approximate transition coupling evaluation.
    [g, h] = approx_tce(P, c, L, T);
    %disp(g(1));
    
    % Entropic transition coupling improvement.
    P = entropic_tci(h, P_old, Px, Py, xi, sink_iter);
end
% In case of numerical instability, make non-negative and normalize.
P = max(P, 0) ./ sum(max(P, 0),2);
if get_sd
    % Get stationary distribution and expected cost of P.
    [stat_dist, exp_cost] = get_best_stat_dist(P, c);
    stat_dist = reshape(stat_dist, dy, dx)';
else
   stat_dist = [];
   exp_cost = g(1);
end
end