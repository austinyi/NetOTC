%% 
% Coordinated optimal transport distance
% 
% Input parameters:
%  -- A1: Symmetric, weighted adjacency matrix for graph 1
%  -- A2: Symmetric, weighted adjacency matrix for graph 2
%  -- sinkiter:  Number of Sinkhorn iterations to perform at each step
%
% Output:
%  -- dist: The computed COPT distance between the two graphs
%  -- P: The COPT plan between the two graphs

function [dist, P] = copt_dist(A1, A2, sinkiter)
    warning('off'); % sqrtm() gives warnings for singular matrices

    m = size(A1, 1);
    n = size(A2, 2);
    L1 = get_lap(A1);
    L2 = get_lap(A2);
    L1inv = pinv(L1);
    L2inv = pinv(L2);
    L2invsqrt = sqrtm(L2inv);

    % Initialize coupling P
    P = normrnd(0, 1, m, n);
    P = m*n*abs(P)./sum(sum(abs(P)));
    logP0 = log(P);
    
    % Function for computing loss
    function loss = get_loss(logP)
        f = zeros(m, 1);
        g = zeros(1, n);
        temp1 = log(n*ones(m, 1));
        temp2 = log(m*ones(1, n));
        for i=1:sinkiter
            %logP = log(n) + logP - (logsumexp(logP, 2));
            %logP = log(m) + logP - (logsumexp(logP, 1));
            
            f = temp1 - logsumexp(logP + g, 2);
            g = temp2 - logsumexp(logP + f, 1);
            logP = f+logP+g;
        end

        sqrt_term = sqrtm(L2invsqrt*exp(logP)'*L1inv*exp(logP)*L2invsqrt);
        loss = m*trace(L1inv) + n*trace(L2inv) - 2*trace(sqrt_term);
        loss = real(loss);
    end

    %options = optimoptions('fminunc','Display','iter-detailed');
    options = optimoptions('fminunc', 'Display', 'off');
    [logP, dist] = fminunc(@(logP)get_loss(logP), logP0, options);
    
    % Compute distance
    dist = sqrt(dist/(m*n));
    
    % Compute optimal P via Sinkhorn iterations
    f = zeros(m, 1);
    g = zeros(1, n);
    for i=1:sinkiter
        f = log(n*ones(m, 1)) - logsumexp(logP + g, 2);
        g = log(m*ones(1, n)) - logsumexp(logP + f, 1);
        logP = f+logP+g;
    end
    P = exp(logP);
end