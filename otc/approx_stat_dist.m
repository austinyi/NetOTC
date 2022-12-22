

function dist = approx_stat_dist(P, iter)
    n = size(P, 1);
    dist = zeros(1, n);
    dist(1) = 1;
    for i=1:iter
        dist = dist*P;
    end
end
