%% 
% K-Nearest Neighbor classification
% 
% Input parameters:
%  -- dist_mat: Distance matrix with rows corresponding to training data
%  and columns corresponding to testing data
%  -- train_classes: Vector of classes for the training data
%  -- k:  Number of neighbors
%
% Output:
%  -- classes: Predicted classes for testing set

function classes = knn(dist_mat, train_labels, k)
    n_test = size(dist_mat, 2);
    classes = zeros(n_test, 1);
    
    for i=1:n_test
        idxs = find(ismember(dist_mat(:,i), min_k(dist_mat(:,i), k)));
        classes(i) = mode(train_labels(idxs));
    end
end