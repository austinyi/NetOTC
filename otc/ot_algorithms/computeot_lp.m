%%
% Computes Optimal Transport using a MATLAB linear programming solver.
%

function [lp_sol,lp_val] = computeot_lp( C,r,c )
% vectorize P and C by: column 1, column 2, etc.
nx = size(r, 1);
ny = size(c, 2);
Aeq = zeros(nx+ny,nx*ny);
beq = [r;c'];

% column sums correct
for row=1:nx
    for t=1:ny
        Aeq(row,(row-1)*ny+t)=1;
    end
end

% row sums correct
for row=nx+1:nx+ny
    for t=0:nx-1
        Aeq(row,(row-nx)+t*ny) = 1;
    end
end

% ensure positivity of each entry
lb = zeros(nx*ny,1);

% solve OT LP using linprog
cost = reshape(C,nx*ny,1);
options = optimoptions('linprog','Display','none');
[lp_sol,lp_val] = linprog(cost,[],[],Aeq,beq,lb,[],options);
end