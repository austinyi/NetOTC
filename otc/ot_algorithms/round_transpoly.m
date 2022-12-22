%%
% Implementation of our algorithm for rounding a matrix onto the U_{r,c}
% transport polytope. See Section 2 of the paper for details.
%

function A = round_transpoly( X,r,c )

A=X;
n1=size(A,1);
n2=size(A,2);
r_A = sum(A,2);
for i=1:n1
    scaling = min(1,r(i)/r_A(i));
    A(i,:)=scaling*A(i,:);
end

c_A = sum(A,1);
for j=1:n2
    scaling = min(1,c(j)/c_A(j));
    A(:,j)=scaling*A(:,j);
end

r_A = sum(A,2);
c_A = sum(A,1);
err_r = r_A - r;
err_c = c_A - c;

if ~all(err_r==0) & ~all(err_c==0)
    A = A + err_r*err_c/sum(abs(err_r));
end
end

