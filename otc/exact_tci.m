%%
% exact_tci.m
%
% Exact transition coupling improvement.

function P = exact_tci(g, h, P0, Px, Py)
x_sizes = size(Px);
y_sizes = size(Py);
dx = x_sizes(1);
dy = y_sizes(1);
P = zeros(dx*dy, dx*dy);
%% Try to improve with respect to g.
% Check if g is constant.
g_const = 1;
for i=1:dx
    for j=(i+1):dy
        if abs(g(i) - g(j)) > 1e-3
            g_const = 0;
            break
        end
    end
end
% If g is not constant, improve transition coupling against g.
if ~g_const
    g_mat = reshape(g, dy, dx)';
    for x_row=1:dx
        for y_row=1:dy
           dist_x = Px(x_row,:);
           dist_y = Py(y_row,:);
           % Check if either distribution is degenerate.
           if any(dist_x == 1) | any(dist_y == 1)
               sol = dist_x' * dist_y;
           % If not degenerate, proceed with OT.
           else
               [sol, val] = computeot_lp(g_mat', dist_x', dist_y);
           end
           idx = dy*(x_row-1)+y_row;
           P(idx,:) = reshape(sol', [], dx*dy);
        end
    end
    if max(abs(P0*g - P*g)) <= 1e-7
        P = P0;
    else
        return
    end
end
    
%% Try to improve with respect to h.
h_mat = reshape(h, dy, dx)';
for x_row=1:dx
    for y_row=1:dy
       dist_x = Px(x_row,:);
       dist_y = Py(y_row,:);
       % Check if either distribution is degenerate.
       if any(dist_x == 1) | any(dist_y == 1)
           sol = dist_x' * dist_y;
       % If not degenerate, proceed with OT.
       else
           [sol, val] = computeot_lp(h_mat', dist_x', dist_y);
       end
       idx = dy*(x_row-1)+y_row;       
       P(idx,:) = reshape(sol', [], dx*dy);
    end
end
if max(abs(P0*h - P*h)) <= 1e-4
    P = P0;
end
end