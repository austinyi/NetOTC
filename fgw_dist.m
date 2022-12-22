%%
%   Fused Gromov-Wasserstein distance.
%

function [FGW,pi] = fgw_dist(M, C1, C2, mu1, mu2, q, alpha)

    % Define some helper functions.
    function loss = fgw_loss(pi)
        loss = (1-alpha)*sum(sum(M.^q.*pi));
        m = size(M,1);
        n = size(M,2);
        for i=1:m
            for j=1:n
                for k=1:m
                    for l=1:n
                        loss = loss + 2*alpha*abs(C1(i,k)-C2(j,l))^q*pi(i,j)*pi(k,l);
                    end
                end
            end
        end
    end

    function grad = fgw_grad(pi)
        grad = (1-alpha)*(M.^q);
        m = size(M,1);
        n = size(M,2);
        for i=1:m
            for j=1:n
                for k=1:m
                    for l=1:n
                        grad(i,j) = grad(i,j) + 2*alpha*abs(C1(i,k)-C2(j,l))^q*pi(k,l);
                    end
                end
            end
        end
    end

    % Initialize coupling
    pi = mu1 .* mu2';
    m = size(pi,1);
    n = size(pi,2);
    
    % Run algorithm
    n_iter = 100;
    for iter=1:n_iter
       % Compute gradient
       G = fgw_grad(pi);
       % Solve OT problem with cost G
       [pi_new, ~] = computeot_lp(G', mu1, mu2');
       pi_new = reshape(pi_new',n,m)';
       % Line search
       fun = @(tau) (fgw_loss((1-tau)*pi+tau*pi_new));
       %tau = fminbnd(fun,0,1);
       tau_vec = 0:0.1:1;
       [~,tau_idx] = min(arrayfun(fun, tau_vec));
       tau = tau_vec(tau_idx);
       % Store updated coupling
       pi = (1-tau)*pi + tau*pi_new;
    end
    
    % Store result
    FGW = fgw_loss(pi);
end