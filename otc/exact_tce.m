%%
% exact_tce.m
%
% Exact transition coupling evaluation.

function [g, h] = exact_tce(P, c)
d = size(P, 1);
c = reshape(c', d, []);
A = [eye(d)-P, zeros(d), zeros(d)
    eye(d), eye(d)-P, zeros(d)
    zeros(d), eye(d), eye(d)-P];
b = [zeros(d, 1); c; zeros(d, 1)];
[sol, r] = linsolve(A, b);
% If the linear solver failed. Use pseudoinverse.
if r == 0
    sol = pinv(A)*b;
end
sol = sol';
g = sol(1:d)';
h = sol((d+1):(2*d))';
end