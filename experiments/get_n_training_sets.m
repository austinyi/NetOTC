%% 
% Get n training / testing sets
% 
% Input parameters:
%  -- size: 
%  -- n:  
%
% Output:
%  -- train_idxs:
%  -- test_idxs: 

function [train_idxs, test_idxs] = get_n_training_sets(size, n)
    p = 0.8; % percent of training data
    
    cv = cvpartition(size,'HoldOut',1-p);
    test_idxs = find(cv.test);
    train_idxs = find(cv.training);
    
    for i=2:n
        cv = cvpartition(size,'HoldOut',1-p);
        test_idxs = cat(2, test_idxs, find(cv.test));
        train_idxs = cat(2, train_idxs, find(cv.training));
    end
end