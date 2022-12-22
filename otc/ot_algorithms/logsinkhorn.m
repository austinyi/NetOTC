%% 
% Implementation of classical Sinkhorn algorithm for matrix scaling.
% Each iteration simply alternately updates (projects) all rows or
% all columns to have correct marginals.
% 
% Input parameters:
%  -- A:  -xi*C
%  -- r:  desired row sums (marginals)         (dims: nx1)
%  -- c:  desired column sums (marginals)      (dims: 1xn)
%  -- T:  number of full Sinkhorn iterations (normalize ALL row or cols)
%  -- C:  cost matrix for OT
%
% Output:
%  -- P:   final scaled matrix

function P = logsinkhorn(A,r,c,T)
[dx, dy] = size(A);
f = zeros(dx, 1);
g = zeros(1, dy);
for t=1:T
    if mod(t,2)==1
        % rescale rows
        f = log(r) - logsumexp(A + g, 2);
    else
        % rescale columns
        g = log(c) - logsumexp(A + f, 1);
    end
end
P = round_transpoly(exp(f+A+g),r,c);
end