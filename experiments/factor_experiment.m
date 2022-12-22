% Network Factors Experiment

% Seed number
rng(311);

% Params
n_iter = 100;
block_size = 5;
n_blocks = 6;
n = block_size*n_blocks;
dim = 5;

for mean1=1:4
    values = table([], []);
    values.Properties.VariableNames = {'Algorithm' 'Accuracy'};
    mean_sigma = 0.5 + mean1 * 0.5;
    disp(["mean_sigma ", num2str(mean_sigma)]);
    for iter=1:n_iter

        %% Setup graph
        % Sample vertices of G1 and G2
        [V1, V2] = get_pointcloud(mean_sigma, block_size, n_blocks, dim);

        A2 = randi(10, n_blocks, n_blocks);

        A1 = zeros(n,n);
        for i = 1:n_blocks
            for j = 1:n_blocks
                for k = 1:block_size
                    rv = rand(1,block_size);
                    rv = rv * A2(i,j) / sum(rv) / block_size;
                    A1((i-1)*block_size+k, (j-1)*block_size+1:j*block_size) = rv;
                end
            end
        end

        %% Do GraphOTC
        % Get transition matrices
        P1 = A1 ./ sum(A1, 2);
        P2 = A2 ./ sum(A2, 2);
        
        stat_dist1 = approx_stat_dist(P1, 100)';
        stat_dist2 = approx_stat_dist(P2, 100)';
        
        stat_dist3 = ones(n,1)/n;
        stat_dist4 = ones(n_blocks,1)/n_blocks;
        
        % Get cost matrix
        c = zeros([n, n_blocks]);
        for i=1:n
            for j=1:n_blocks
                c(i, j) = sum((V1(i,:)-V2(j,:)).^2);
            end
        end

        % Run algorithms
        [~, otc_edge_alignment, otc_alignment] = exact_otc(P1, P2, c);
        [~, fgw_alignment] = fgw_dist(c, A1, A2, stat_dist3, stat_dist4, 1, 0.5);
        [otsd_alignment, ~] = computeot_lp(c', stat_dist1, stat_dist2');
        otsd_alignment = reshape(otsd_alignment, n_blocks, n)';

        aligned_mass_otc = eval_alignment(otc_alignment, block_size, n_blocks);
        aligned_mass_fgw = eval_alignment(fgw_alignment, block_size, n_blocks);
        aligned_mass_otsd = eval_alignment(otsd_alignment, block_size, n_blocks);

        % Store results
        values = [values; {'OTC' aligned_mass_otc}];
        values = [values; {'FGW' aligned_mass_fgw}];          
        values = [values; {'OT-SD' aligned_mass_otsd}];
    end

    disp(['OTC mean accuracy: ' num2str(mean(values{strcmp(values{:,1}, 'OTC'),2}))]);
    disp(['FGW mean accuracy: ' num2str(mean(values{strcmp(values{:,1}, 'FGW'),2}))]);
    disp(['OT-SD mean accuracy: ' num2str(mean(values{strcmp(values{:,1}, 'OT-SD'),2}))]);

    disp(['OTC accuracy standard deviation: ' num2str(std(values{strcmp(values{:,1}, 'OTC'),2}))]);
    disp(['FGW accuracy standard deviation: ' num2str(std(values{strcmp(values{:,1}, 'FGW'),2}))]);
    disp(['OT-SD accuracy standard deviation: ' num2str(std(values{strcmp(values{:,1}, 'OT-SD'),2}))]);
end



%% Helper functions
% Point cloud sampling function
function [points, factor] = get_pointcloud(mean_sigma, n_points, n_clouds, dimension)
    point_sigma = 1;
    
    % Randomly draw means of each Gaussian
    means = normrnd(0, mean_sigma, n_clouds, dimension);
    factor = means;
    
    % Draw point clouds
    points = zeros(n_points*n_clouds, dimension);
    for c = 1:n_clouds
        for dim=1:dimension
            points(((c-1)*n_points+1):c*n_points, dim) = normrnd(means(c,dim), point_sigma, [n_points, 1]);
        end
    end
end

% Function for adding up mass in correct alignment
function alignment = eval_alignment(coupling, block_size, n_blocks)
    alignment = 0;
    for b=1:n_blocks
        for i=1:block_size
            idx = (b-1)*block_size + i;
            alignment = alignment + coupling(idx, b);
        end
    end
end
