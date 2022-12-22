%%
% pc_align_cv.m
%%

function fgw_params = sbm_align_cv(n_points1, n_points2, cv_iter, alpha_grid)

    n_clouds = 4;

    % Setup variables
    fgw_acc = zeros(length(alpha_grid), cv_iter);
    fgw_edge_acc = zeros(length(alpha_grid), cv_iter);

    %parpool(5);    
    %parfor iter=1:cv_iter
    for iter=1:cv_iter
        disp(['CV Iteration' num2str(iter)]);
        
        % Iteration accuracy
        iter_fgw_acc = zeros(length(alpha_grid), 1);
        iter_fgw_edge_acc = zeros(length(alpha_grid), 1);
        

        % Generate SBM
        [~,A1] = stochastic_block_model([n_points1 n_points1 n_points1 n_points1], [1.0 0.1 0.1 0.1; 0.1 0.8 0.1 0.1; 0.1 0.1 0.6 0.1; 0.1 0.1 0.1 0.4]);
        [~,A2] = stochastic_block_model([n_points2 n_points2 n_points2 n_points2], [1.0 0.1 0.1 0.1; 0.1 0.8 0.1 0.1; 0.1 0.1 0.6 0.1; 0.1 0.1 0.1 0.4]);

        n1 = size(A1,1);
        n2 = size(A2,1);

        unif_dist1 = ones(n1,1)/n1;
        unif_dist2 = ones(n2,1)/n2;

        % Construct cost function
        c = standardized_degree_cost(A1, A2);


        for alpha_idx = 1:length(alpha_grid)
            alpha = alpha_grid(alpha_idx);
            
            % Run FGW and evaluate
            [~, fgw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, alpha);

            aligned_mass_fgw = 0;
            for cl = 1:n_clouds
                aligned_mass_fgw = aligned_mass_fgw + sum(fgw_alignment(((cl-1)*n_points1 + 1):cl*n_points1, ((cl-1)*n_points2 + 1):cl*n_points2), 'all');
            end            
            iter_fgw_acc(alpha_idx) = aligned_mass_fgw;
            disp(['FGW Accuracy: ' num2str(aligned_mass_fgw)]);

            edge_mass_fgw = 0;
            for x_row=1:n1
                for x_col=1:n1
                    for y_row=1:n2
                        for y_col=1:n2
                            if fix((x_row-1)/n_points1) == fix((x_col-1)/n_points1) && fix((y_row-1)/n_points2) == fix((y_col-1)/n_points2)
                                edge_mass_fgw = edge_mass_fgw + fgw_alignment(x_row, y_row)*fgw_alignment(x_col, y_col);
                            end
                        end
                    end
                end
            end
            iter_fgw_edge_acc(alpha_idx) = edge_mass_fgw;
            disp(['FGW Edge Accuracy: ' num2str(edge_mass_fgw)]);
        end
        
        fgw_acc(:,iter) = iter_fgw_acc;
        fgw_edge_acc(:,iter) = iter_fgw_edge_acc;
    end
    %delete(gcp('nocreate'));
    
    fgw_acc = mean(fgw_acc, 2);
    disp('FGW mean accuracy');
    disp(fgw_acc);
    
    fgw_edge_acc = mean(fgw_edge_acc, 2);
    disp('FGW mean edge accuracy');
    disp(fgw_edge_acc);
  
    alpha_idx1 = find(fgw_acc == max(max(fgw_acc)));
    alpha_idx2 = find(fgw_edge_acc == max(max(fgw_edge_acc)));
    fgw_params = [alpha_grid(alpha_idx1(1)); alpha_grid(alpha_idx2(1))];
end