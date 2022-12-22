%%
% sbm_alignment_experiment.m
%%

function sbm_alignment(n_points1, n_points2, n_iter, longleaf)
    rng(1000);
    
    % Setup directories
    expid = strrep(num2str(now, '%f'), '.', '_');
    disp(expid)
    
    format longG;
    if longleaf
        addpath('/nas/longleaf/home/bongsoo/GraphOTC_final/GraphOTC/','/nas/longleaf/home/bongsoo/GraphOTC_final/GraphOTC/otc/','/nas/longleaf/home/bongsoo/GraphOTC_final/GraphOTC/otc/ot_algorithms/','/nas/longleaf/home/bongsoo/GraphOTC_final/GraphOTC/utils/', '/nas/longleaf/home/bongsoo/GraphOTC_final/GraphOTC/experiments/');
        dir = '/pine/scr/b/o/bongsoo/GraphOTC/sbm_alignment/';
    else
        dir = '/Users/bongsooyi/Documents/NetOTC/experiments/isomorphism/';
    end
    
    expdir = [dir expid '/'];
    mkdir(expdir);

   

    % Experiment parameters
    n_clouds = 4;


    % Cross-validation parameters
    do_cv = 1;
    cv_iter = 10;
    alpha_grid = 0:0.1:1;

    % Setup data tables
    values = table([], []);
    values.Properties.VariableNames = {'Algorithm' 'Accuracy'};
    edge_values = table([], []);
    edge_values.Properties.VariableNames = {'Algorithm' 'Accuracy'};



    
    if do_cv
        disp('Doing cross-validation');
        fgw_params = sbm_align_cv(n_points1, n_points2, cv_iter, alpha_grid);
        disp('Best FGW params:');
        disp(fgw_params);
    else
        fgw_params = [0.5, 0.5];
    end


    iter = 0;

    while iter < n_iter
        try 
            % Generate SBM
            [~,A1] = stochastic_block_model([n_points1 n_points1 n_points1 n_points1], [1.0 0.1 0.1 0.1; 0.1 0.8 0.1 0.1; 0.1 0.1 0.6 0.1; 0.1 0.1 0.1 0.4]);
            [~,A2] = stochastic_block_model([n_points2 n_points2 n_points2 n_points2], [1.0 0.1 0.1 0.1; 0.1 0.8 0.1 0.1; 0.1 0.1 0.6 0.1; 0.1 0.1 0.1 0.4]);
            
            n1 = size(A1,1);
            n2 = size(A2,1);

            % Compute transition matrices and stationary distributions
            P1 = A1 ./ sum(A1, 2);
            P2 = A2 ./ sum(A2, 2);

            stat_dist1 = approx_stat_dist(P1, 100)';
            stat_dist2 = approx_stat_dist(P2, 100)';

            unif_dist1 = ones(n1,1)/n1;
            unif_dist2 = ones(n2,1)/n2;
          
            % Construct cost function
            %c = get_degree_cost(A1, A2);
            c = standardized_degree_cost(A1, A2);
            
            % Run algorithms
            [~, otc_edge_alignment, otc_alignment] = exact_otc(P1, P2, c);
            [~, fgw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, fgw_params(1));
            [~, gw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 1);
            [otsd_alignment, ~] = computeot_lp(c', stat_dist1, stat_dist2');
            otsd_alignment = reshape(otsd_alignment, n2, n1)';
            [~, copt_alignment] = copt_dist(A1, A2, 10);
            copt_alignment = copt_alignment ./ sum(copt_alignment, 'all');


            iter = iter + 1;
            disp(['Iteration' num2str(iter)]);          

            %% Evaluate the alignments
            % Vertex alignment
            aligned_mass_otc = eval_alignment(otc_alignment);
            aligned_mass_fgw = eval_alignment(fgw_alignment);
            aligned_mass_gw = eval_alignment(gw_alignment);
            aligned_mass_otsd = eval_alignment(otsd_alignment);
            aligned_mass_copt = eval_alignment(copt_alignment);

            disp(['OTC Accuracy: ' num2str(aligned_mass_otc)]);
            disp(['FGW Accuracy: ' num2str(aligned_mass_fgw)]);
            disp(['GW Accuracy: ' num2str(aligned_mass_gw)]);            
            disp(['OT-SD Accuracy: ' num2str(aligned_mass_otsd)]);
            disp(['COPT Accuracy: ' num2str(aligned_mass_copt)]);

            values = [values; {'OTC' aligned_mass_otc}];
            values = [values; {'FGW' aligned_mass_fgw}];
            values = [values; {'GW' aligned_mass_gw}];            
            values = [values; {'OT-SD' aligned_mass_otsd}];
            values = [values; {'COPT' aligned_mass_copt}];

            % Edge alignment
            % Refit FGW with edge-tuned parameters
            [~, fgw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, fgw_params(2));            

            edge_mass_otc = 0;
            for x_row=1:n1
                for x_col=1:n1
                    for y_row=1:n2
                        for y_col=1:n2
                            idx1 = n2*(x_row-1)+y_row;
                            idx2 = n2*(x_col-1)+y_col;
                            if fix((x_row-1)/n_points1) == fix((x_col-1)/n_points1) && fix((y_row-1)/n_points2) == fix((y_col-1)/n_points2)
                                edge_mass_otc = edge_mass_otc + otc_edge_alignment(idx1, idx2)*otc_alignment(x_row, y_row);
                            end
                        end
                    end
                end
            end

            edge_mass_fgw = eval_edge_alignment(fgw_alignment);
            edge_mass_gw = eval_edge_alignment(gw_alignment);
            edge_mass_otsd = eval_edge_alignment(otsd_alignment);
            edge_mass_copt = eval_edge_alignment(copt_alignment);

            disp(['OTC Edge Accuracy: ' num2str(edge_mass_otc)]);
            disp(['FGW Edge Accuracy: ' num2str(edge_mass_fgw)]);
            disp(['GW Edge Accuracy: ' num2str(edge_mass_gw)]);     
            disp(['OT-SD Edge Accuracy: ' num2str(edge_mass_otsd)]); 
            disp(['COPT Edge Accuracy: ' num2str(edge_mass_copt)]); 

            edge_values = [edge_values; {'OTC' edge_mass_otc}];
            edge_values = [edge_values; {'FGW' edge_mass_fgw}];
            edge_values = [edge_values; {'GW' edge_mass_gw}];
            edge_values = [edge_values; {'OT-SD' edge_mass_otsd}];
            edge_values = [edge_values; {'COPT' edge_mass_copt}];

        catch
        end
    end
    %delete(gcp('nocreate'));


    % Write accuracies
    %disp(values);
    writetable(values, [expdir 'alignment_accuracies.txt']);
    %disp(edge_values);
    writetable(edge_values, [expdir 'edge_alignment_accuracies.txt']);

    % Print mean accuracy 
    disp(['OTC mean accuracy: ' num2str(mean(values{strcmp(values{:,1}, 'OTC'),2}))]);
    disp(['FGW mean accuracy: ' num2str(mean(values{strcmp(values{:,1}, 'FGW'),2}))]);
    disp(['GW mean accuracy: ' num2str(mean(values{strcmp(values{:,1}, 'GW'),2}))]);
    disp(['OT-SD mean accuracy: ' num2str(mean(values{strcmp(values{:,1}, 'OT-SD'),2}))]);
    disp(['COPT mean accuracy: ' num2str(mean(values{strcmp(values{:,1}, 'COPT'),2}))]);

    disp(['OTC edge mean accuracy: ' num2str(mean(edge_values{strcmp(edge_values{:,1}, 'OTC'),2}))]);
    disp(['FGW edge mean accuracy: ' num2str(mean(edge_values{strcmp(edge_values{:,1}, 'FGW'),2}))]);
    disp(['GW edge mean accuracy: ' num2str(mean(edge_values{strcmp(edge_values{:,1}, 'GW'),2}))]);
    disp(['OT-SD edge mean accuracy: ' num2str(mean(edge_values{strcmp(edge_values{:,1}, 'OT-SD'),2}))]);
    disp(['COPT edge mean accuracy: ' num2str(mean(edge_values{strcmp(edge_values{:,1}, 'COPT'),2}))]);

    % Vertex alignment evaluation function
    function score = eval_alignment(alignment)
        %global n_clouds n_points1 n_points2;
        score = 0;
        for cl = 1:n_clouds
            score = score + sum(alignment(((cl-1)*n_points1 + 1):cl*n_points1, ((cl-1)*n_points2 + 1):cl*n_points2), 'all');
        end            
    end

    % Edge alignment evaluation function for FGW, GW, OT-SD, COPT
    function edge_score = eval_edge_alignment(alignment)
        %global n_points1 n_points2;
        edge_score = 0;
        n1 = size(alignment, 1);
        n2 = size(alignment, 2);
        for x_row=1:n1
            for x_col=1:n1
                for y_row=1:n2
                    for y_col=1:n2
                        if fix((x_row-1)/n_points1) == fix((x_col-1)/n_points1) && fix((y_row-1)/n_points2) == fix((y_col-1)/n_points2)
                            edge_score = edge_score + alignment(x_row, y_row)*alignment(x_col, y_col);
                        end
                    end
                end
            end
        end
    end
end
