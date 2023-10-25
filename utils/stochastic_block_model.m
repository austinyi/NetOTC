%%
% stochastic_block_model.m
%
% Input: sizes, probs
% Output: Graph G, Adjacency matrix A

function [G, A] = stochastic_block_model(sizes, probs)
    
    % Check input type
    if isvector(sizes) == 0
        error('''sizes'' must be a vector.')
    elseif ~ismatrix(probs) || (size(probs, 1) ~= size(probs, 2))
        error('''probs'' must be a sqaure matrix.')
    elseif ~issymmetric(probs)
        error('''probs'' must be symmetric.')
    elseif length(sizes) ~= length(probs(1,:))
        error('''sizes'' and ''p'' do not match.')
    end 
    
   
    n = sum(sizes); % Total number of nodes
    n_b = length(sizes); % Total number of blocks
    A = zeros(n); 
    
    % Column index of each blocks starting
    start = [1];
    s = 1;
    for i = 1:n_b
        s = s + sizes(i);
        start = [start, s];
    end
    
    % Generating Adjacency Matrix (upper)
    %
    % First, generate blocks on diagonal
    for i = 1:n_b
        p = probs(i,i);
        for j = start(i):start(i+1)-1
            for k = j+1:start(i+1)-1
                A(j,k) = randsrc(1, 1, [0, 1; 1 - p, p]);  
            end
        end
    end
    
    % Nondiagonal blocks
    for i = 1:n_b-1
        for j = i+1:n_b
            A(start(i):start(i+1)-1, start(j):start(j+1)-1) = randsrc(sizes(i), sizes(j), [0, 1; 1 - probs(i,j), probs(i,j)]);
        end
    end
    
    % Fill lower traingular matrix
    A = A + A';
    
    % Graph
    G = graph(A);
    
end


 
    