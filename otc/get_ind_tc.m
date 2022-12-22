%%
% get_ind_tc.m
%
% Compute independent transition coupling of two transition matrices.

function [P_ind] = get_ind_tc(Px, Py)
    [dx, dx_col] = size(Px);
    [dy, dy_col] = size(Py);
    
    P_ind = zeros(dx*dy, dx_col*dy_col);
    for x_row=1:dx
        for x_col=1:dx_col
            for y_row=1:dy
                for y_col=1:dy_col
                    idx1 = dy*(x_row-1)+y_row;
                    idx2 = dy*(x_col-1)+y_col;
                    P_ind(idx1, idx2) = Px(x_row, x_col)*Py(y_row, y_col);
                end
            end
        end
    end
end