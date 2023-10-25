function iso = check_isomorphism(idx, A1, A2)
    if size(A1,1) ~= size(A2,1)
        iso = 0;
        return;
    end
    
    if isequal(unique(idx).', 1:size(A1,1)) == 0
        iso = 0;
        return;
    end
    
    [row,col] = find(tril(A1)~=0);
    n = size(row);
    for i=1:n
        if A2(idx(row(i)),idx(col(i))) ~= A1(row(i),col(i))
            iso = 0;
            return;
        end
    end
    
    [row,col] = find(tril(A2)~=0);
    inv_idx = zeros(size(A1,1),1);
    for i=1:size(A1,1)
        inv_idx(i) = find(idx==i);
    end
    n = size(row);
    for i=1:n
        if A1(inv_idx(row(i)),inv_idx(col(i))) ~= A2(row(i),col(i))
            iso = 0;
            return;
        end
    end
    
    iso = 1;
    
end