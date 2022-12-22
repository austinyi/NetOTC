function isomorphism_experiment(n_iter, longleaf)
    rng(1);
    
    % Setup directories
    expid = strrep(num2str(now, '%f'), '.', '_');
    disp(expid)
    
    format longG;
    if longleaf
        addpath('/nas/longleaf/home/bongsoo/GraphOTC_final/GraphOTC/','/nas/longleaf/home/bongsoo/GraphOTC_final/GraphOTC/otc/','/nas/longleaf/home/bongsoo/GraphOTC_final/GraphOTC/otc/ot_algorithms/','/nas/longleaf/home/bongsoo/GraphOTC_final/GraphOTC/utils/', '/nas/longleaf/home/bongsoo/GraphOTC_final/GraphOTC/experiments/',  '/nas/longleaf/home/bongsoo/GraphOTC_final/GraphOTC/examples/');
        dir = '/pine/scr/b/o/bongsoo/GraphOTC/isomorphism/';
    else
        dir = '/Users/bongsooyi/Documents/NetOTC/experiments/isomorphism/';
    end
    
    expdir = [dir expid '/'];
    mkdir(expdir);
    
  
    % 1.
    % Erdos-Renyi random graph
    disp('Erdos-Renyi random graph start (n=6~15, p=1/3)')

    % Number of success for each algorithm
    otc_success = 0;
    otsd_success = 0;
    gw_success = 0;
    fgw_success = 0;
    
    % Number of connected graphs among the randomly generated graphs
    connected = 0;
    
    parpool(24);
    parfor iter=1:n_iter
        try 
            n = 5 + randi(10);
            p = 1/3;
            
            % Construct adjacency matrices
            [~,A1] = stochastic_block_model([n], [p]);
            
            % Random permutation
            perm = randperm(n);
            A2 = A1(perm,perm);

            % Get transition matrices
            P1 = A1 ./ sum(A1, 2);
            P2 = A2 ./ sum(A2, 2);

            % Get distributions
            stat_dist1 = approx_stat_dist(P1, 100)';
            stat_dist2 = approx_stat_dist(P2, 100)';
            unif_dist1 = ones(n,1)/n;
            unif_dist2 = ones(n,1)/n;

            % Get cost function
            c = get_degree_cost(A1, A2);
    
            % Run algorithm
            [~, ~, otc_alignment] = exact_otc(P1, P2, c);
            [otsd_alignment, ~] = computeot_lp(c', stat_dist1, stat_dist2');
            otsd_alignment = reshape(otsd_alignment, n, n)';
            [~, gw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 1);
            [~, fgw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 0.5);

            % Get alignment
            [~,idx_otc] = max(otc_alignment,[],2);
            [~,idx_otsd] = max(otsd_alignment,[],2);
            [~,idx_gw] = max(gw_alignment,[],2);
            [~,idx_fgw] = max(fgw_alignment,[],2);

            % Check success
            otc_success = otc_success + check_isomorphism(idx_otc, A1, A2);
            otsd_success = otsd_success + check_isomorphism(idx_otsd, A1, A2);
            gw_success = gw_success + check_isomorphism(idx_gw, A1, A2);
            fgw_success = fgw_success + check_isomorphism(idx_fgw, A1, A2);
            
            % Update connected graph number
            connected = connected + 1;
        catch
        end
    end
    delete(gcp('nocreate'));
    
    % Display accuracy
    disp(['OTC mean accuracy: ', num2str(otc_success/connected)]);
    disp(['OTSD mean accuracy: ', num2str(otsd_success/connected)]);
    disp(['GW mean accuracy: ', num2str(gw_success/connected)]);
    disp(['FGW mean accuracy: ', num2str(fgw_success/connected)]);

    disp(connected);
  
    
    % 2.
    % Erdos-Renyi random graph
    disp('Erdos-Renyi random graph start (n=6~15, p=2/3)')

    % Number of success for each algorithm
    otc_success = 0;
    otsd_success = 0;
    gw_success = 0;
    fgw_success = 0;
    
    % Number of connected graphs among the randomly generated graphs
    connected = 0;
    
    parpool(24);
    parfor iter=1:n_iter
        try 
            n = 5 + randi(10);
            p = 2/3;
            
            % Construct adjacency matrices
            [~,A1] = stochastic_block_model([n], [p]);
            
            % Random permutation
            perm = randperm(n);
            A2 = A1(perm,perm);

            % Get transition matrices
            P1 = A1 ./ sum(A1, 2);
            P2 = A2 ./ sum(A2, 2);

            % Get distributions
            stat_dist1 = approx_stat_dist(P1, 100)';
            stat_dist2 = approx_stat_dist(P2, 100)';
            unif_dist1 = ones(n,1)/n;
            unif_dist2 = ones(n,1)/n;

            % Get cost function
            c = get_degree_cost(A1, A2);
    
            % Run algorithm
            [~, ~, otc_alignment] = exact_otc(P1, P2, c);
            [otsd_alignment, ~] = computeot_lp(c', stat_dist1, stat_dist2');
            otsd_alignment = reshape(otsd_alignment, n, n)';
            [~, gw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 1);
            [~, fgw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 0.5);

            % Get alignment
            [~,idx_otc] = max(otc_alignment,[],2);
            [~,idx_otsd] = max(otsd_alignment,[],2);
            [~,idx_gw] = max(gw_alignment,[],2);
            [~,idx_fgw] = max(fgw_alignment,[],2);

            % Check success
            otc_success = otc_success + check_isomorphism(idx_otc, A1, A2);
            otsd_success = otsd_success + check_isomorphism(idx_otsd, A1, A2);
            gw_success = gw_success + check_isomorphism(idx_gw, A1, A2);
            fgw_success = fgw_success + check_isomorphism(idx_fgw, A1, A2);
            
            % Update connected graph number
            connected = connected + 1;
        catch
        end
    end
    delete(gcp('nocreate'));
    
    % Display accuracy
    disp(['OTC mean accuracy: ', num2str(otc_success/connected)]);
    disp(['OTSD mean accuracy: ', num2str(otsd_success/connected)]);
    disp(['GW mean accuracy: ', num2str(gw_success/connected)]);
    disp(['FGW mean accuracy: ', num2str(fgw_success/connected)]);

    disp(connected);
    
    
    % 3.
    % Erdos-Renyi random graph
    disp('Erdos-Renyi random graph start (n=16~25, p=1/4)')

    % Number of success for each algorithm
    otc_success = 0;
    otsd_success = 0;
    gw_success = 0;
    fgw_success = 0;
    
    % Number of connected graphs among the randomly generated graphs
    connected = 0;
    
    parpool(24);
    parfor iter=1:n_iter
        try 
            n = 15 + randi(10);
            p = 1/4;
            
            % Construct adjacency matrices
            [~,A1] = stochastic_block_model([n], [p]);
            
            % Random permutation
            perm = randperm(n);
            A2 = A1(perm,perm);

            % Get transition matrices
            P1 = A1 ./ sum(A1, 2);
            P2 = A2 ./ sum(A2, 2);

            % Get distributions
            stat_dist1 = approx_stat_dist(P1, 100)';
            stat_dist2 = approx_stat_dist(P2, 100)';
            unif_dist1 = ones(n,1)/n;
            unif_dist2 = ones(n,1)/n;

            % Get cost function
            c = get_degree_cost(A1, A2);
    
            % Run algorithm
            [~, ~, otc_alignment] = exact_otc(P1, P2, c);
            [otsd_alignment, ~] = computeot_lp(c', stat_dist1, stat_dist2');
            otsd_alignment = reshape(otsd_alignment, n, n)';
            [~, gw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 1);
            [~, fgw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 0.5);

            % Get alignment
            [~,idx_otc] = max(otc_alignment,[],2);
            [~,idx_otsd] = max(otsd_alignment,[],2);
            [~,idx_gw] = max(gw_alignment,[],2);
            [~,idx_fgw] = max(fgw_alignment,[],2);

            % Check success
            otc_success = otc_success + check_isomorphism(idx_otc, A1, A2);
            otsd_success = otsd_success + check_isomorphism(idx_otsd, A1, A2);
            gw_success = gw_success + check_isomorphism(idx_gw, A1, A2);
            fgw_success = fgw_success + check_isomorphism(idx_fgw, A1, A2);
            
            % Update connected graph number
            connected = connected + 1;
        catch
        end
    end
    delete(gcp('nocreate'));
    
    % Display accuracy
    disp(['OTC mean accuracy: ', num2str(otc_success/connected)]);
    disp(['OTSD mean accuracy: ', num2str(otsd_success/connected)]);
    disp(['GW mean accuracy: ', num2str(gw_success/connected)]);
    disp(['FGW mean accuracy: ', num2str(fgw_success/connected)]);

    disp(connected);
  
    
    % 4.
    % Erdos-Renyi random graph
    disp('Erdos-Renyi random graph start (n=16~25, p=3/4)')

    % Number of success for each algorithm
    otc_success = 0;
    otsd_success = 0;
    gw_success = 0;
    fgw_success = 0;
    
    % Number of connected graphs among the randomly generated graphs
    connected = 0;
    
    parpool(24);
    parfor iter=1:n_iter
        try 
            n = 15 + randi(10);
            p = 3/4;
            
            % Construct adjacency matrices
            [~,A1] = stochastic_block_model([n], [p]);
            
            % Random permutation
            perm = randperm(n);
            A2 = A1(perm,perm);

            % Get transition matrices
            P1 = A1 ./ sum(A1, 2);
            P2 = A2 ./ sum(A2, 2);

            % Get distributions
            stat_dist1 = approx_stat_dist(P1, 100)';
            stat_dist2 = approx_stat_dist(P2, 100)';
            unif_dist1 = ones(n,1)/n;
            unif_dist2 = ones(n,1)/n;

            % Get cost function
            c = get_degree_cost(A1, A2);
    
            % Run algorithm
            [~, ~, otc_alignment] = exact_otc(P1, P2, c);
            [otsd_alignment, ~] = computeot_lp(c', stat_dist1, stat_dist2');
            otsd_alignment = reshape(otsd_alignment, n, n)';
            [~, gw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 1);
            [~, fgw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 0.5);

            % Get alignment
            [~,idx_otc] = max(otc_alignment,[],2);
            [~,idx_otsd] = max(otsd_alignment,[],2);
            [~,idx_gw] = max(gw_alignment,[],2);
            [~,idx_fgw] = max(fgw_alignment,[],2);

            % Check success
            otc_success = otc_success + check_isomorphism(idx_otc, A1, A2);
            otsd_success = otsd_success + check_isomorphism(idx_otsd, A1, A2);
            gw_success = gw_success + check_isomorphism(idx_gw, A1, A2);
            fgw_success = fgw_success + check_isomorphism(idx_fgw, A1, A2);
            
            % Update connected graph number
            connected = connected + 1;
        catch
        end
    end
    delete(gcp('nocreate'));
    
    % Display accuracy
    disp(['OTC mean accuracy: ', num2str(otc_success/connected)]);
    disp(['OTSD mean accuracy: ', num2str(otsd_success/connected)]);
    disp(['GW mean accuracy: ', num2str(gw_success/connected)]);
    disp(['FGW mean accuracy: ', num2str(fgw_success/connected)]);

    disp(connected);


    % 5.
    % SBM 7-7-7-7
    disp('SBM 7-7-7-7 start')

    % Number of success for each algorithm
    otc_success = 0;
    otsd_success = 0;
    gw_success = 0;
    fgw_success = 0;
    
    % Number of connected graphs among the randomly generated graphs
    connected = 0;
    
    parpool(24);
    parfor iter=1:n_iter
        try 
            % Construct adjacency matrices
            [~,A1] = stochastic_block_model([7 7 7 7], [0.7 0.1 0.1 0.1; 0.1 0.7 0.1 0.1; 0.1 0.1 0.7 0.1; 0.1 0.1 0.1 0.7]);
            n = length(A1);
            
            % Random permutation
            perm = randperm(n);
            A2 = A1(perm,perm);

            % Get transition matrices
            P1 = A1 ./ sum(A1, 2);
            P2 = A2 ./ sum(A2, 2);

            % Get distributions
            stat_dist1 = approx_stat_dist(P1, 100)';
            stat_dist2 = approx_stat_dist(P2, 100)';
            unif_dist1 = ones(n,1)/n;
            unif_dist2 = ones(n,1)/n;

            % Get cost function
            c = get_degree_cost(A1, A2);
    
            % Run algorithm
            [~, ~, otc_alignment] = exact_otc(P1, P2, c);
            [otsd_alignment, ~] = computeot_lp(c', stat_dist1, stat_dist2');
            otsd_alignment = reshape(otsd_alignment, n, n)';
            [~, gw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 1);
            [~, fgw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 0.5);

            % Get alignment
            [~,idx_otc] = max(otc_alignment,[],2);
            [~,idx_otsd] = max(otsd_alignment,[],2);
            [~,idx_gw] = max(gw_alignment,[],2);
            [~,idx_fgw] = max(fgw_alignment,[],2);

            % Check success
            otc_success = otc_success + check_isomorphism(idx_otc, A1, A2);
            otsd_success = otsd_success + check_isomorphism(idx_otsd, A1, A2);
            gw_success = gw_success + check_isomorphism(idx_gw, A1, A2);
            fgw_success = fgw_success + check_isomorphism(idx_fgw, A1, A2);
            
            % Update connected graph number
            connected = connected + 1;
        catch
        end
    end
    delete(gcp('nocreate'));
    
    % Display accuracy
    disp(['OTC mean accuracy: ', num2str(otc_success/connected)]);
    disp(['OTSD mean accuracy: ', num2str(otsd_success/connected)]);
    disp(['GW mean accuracy: ', num2str(gw_success/connected)]);
    disp(['FGW mean accuracy: ', num2str(fgw_success/connected)]);

    disp(connected);


    % 6.
    % SBM 10-8-6
    disp('SBM 10-8-6 start')

    % Number of success for each algorithm
    otc_success = 0;
    otsd_success = 0;
    gw_success = 0;
    fgw_success = 0;
    
    % Number of connected graphs among the randomly generated graphs
    connected = 0;
    
    parpool(24);
    parfor iter=1:n_iter
        try 
            % Construct adjacency matrices
            [~,A1] = stochastic_block_model([10 8 6], [0.7 0.1 0.1; 0.1 0.7 0.1; 0.1 0.1 0.7]);
            n = length(A1);
            
            % Random permutation
            perm = randperm(n);
            A2 = A1(perm,perm);

            % Get transition matrices
            P1 = A1 ./ sum(A1, 2);
            P2 = A2 ./ sum(A2, 2);

            % Get distributions
            stat_dist1 = approx_stat_dist(P1, 100)';
            stat_dist2 = approx_stat_dist(P2, 100)';
            unif_dist1 = ones(n,1)/n;
            unif_dist2 = ones(n,1)/n;

            % Get cost function
            c = get_degree_cost(A1, A2);
    
            % Run algorithm
            [~, ~, otc_alignment] = exact_otc(P1, P2, c);
            [otsd_alignment, ~] = computeot_lp(c', stat_dist1, stat_dist2');
            otsd_alignment = reshape(otsd_alignment, n, n)';
            [~, gw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 1);
            [~, fgw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 0.5);

            % Get alignment
            [~,idx_otc] = max(otc_alignment,[],2);
            [~,idx_otsd] = max(otsd_alignment,[],2);
            [~,idx_gw] = max(gw_alignment,[],2);
            [~,idx_fgw] = max(fgw_alignment,[],2);

            % Check success
            otc_success = otc_success + check_isomorphism(idx_otc, A1, A2);
            otsd_success = otsd_success + check_isomorphism(idx_otsd, A1, A2);
            gw_success = gw_success + check_isomorphism(idx_gw, A1, A2);
            fgw_success = fgw_success + check_isomorphism(idx_fgw, A1, A2);
            
            % Update connected graph number
            connected = connected + 1;
        catch
        end
    end
    delete(gcp('nocreate'));
    
    % Display accuracy
    disp(['OTC mean accuracy: ', num2str(otc_success/connected)]);
    disp(['OTSD mean accuracy: ', num2str(otsd_success/connected)]);
    disp(['GW mean accuracy: ', num2str(gw_success/connected)]);
    disp(['FGW mean accuracy: ', num2str(fgw_success/connected)]);

    disp(connected);
    
    
    % 7.
    % SBM 7-7-7
    disp('SBM 7-7-7 start')

    % Number of success for each algorithm
    otc_success = 0;
    otsd_success = 0;
    gw_success = 0;
    fgw_success = 0;
    
    % Number of connected graphs among the randomly generated graphs
    connected = 0;
    
    parpool(24);
    parfor iter=1:n_iter
        try 
            % Construct adjacency matrices
            [~,A1] = stochastic_block_model([7 7 7], [0.7 0.1 0.1; 0.1 0.7 0.1; 0.1 0.1 0.7]);
            n = length(A1);
            
            % Random permutation
            perm = randperm(n);
            A2 = A1(perm,perm);

            % Get transition matrices
            P1 = A1 ./ sum(A1, 2);
            P2 = A2 ./ sum(A2, 2);

            % Get distributions
            stat_dist1 = approx_stat_dist(P1, 100)';
            stat_dist2 = approx_stat_dist(P2, 100)';
            unif_dist1 = ones(n,1)/n;
            unif_dist2 = ones(n,1)/n;

            % Get cost function
            c = get_degree_cost(A1, A2);
    
            % Run algorithm
            [~, ~, otc_alignment] = exact_otc(P1, P2, c);
            [otsd_alignment, ~] = computeot_lp(c', stat_dist1, stat_dist2');
            otsd_alignment = reshape(otsd_alignment, n, n)';
            [~, gw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 1);
            [~, fgw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 0.5);

            % Get alignment
            [~,idx_otc] = max(otc_alignment,[],2);
            [~,idx_otsd] = max(otsd_alignment,[],2);
            [~,idx_gw] = max(gw_alignment,[],2);
            [~,idx_fgw] = max(fgw_alignment,[],2);

            % Check success
            otc_success = otc_success + check_isomorphism(idx_otc, A1, A2);
            otsd_success = otsd_success + check_isomorphism(idx_otsd, A1, A2);
            gw_success = gw_success + check_isomorphism(idx_gw, A1, A2);
            fgw_success = fgw_success + check_isomorphism(idx_fgw, A1, A2);
            
            % Update connected graph number
            connected = connected + 1;
        catch
        end
    end
    delete(gcp('nocreate'));
    
    % Display accuracy
    disp(['OTC mean accuracy: ', num2str(otc_success/connected)]);
    disp(['OTSD mean accuracy: ', num2str(otsd_success/connected)]);
    disp(['GW mean accuracy: ', num2str(gw_success/connected)]);
    disp(['FGW mean accuracy: ', num2str(fgw_success/connected)]);

    disp(connected);
    
    
    % 8.
    % Random weighted adjacency matrix (0,1,2)
    disp('Random weighted adjacency matrix (0,1,2) start')
    
    % Number of success for each algorithm
    otc_success = 0;
    otsd_success = 0;
    gw_success = 0;
    fgw_success = 0;
    
    % Number of connected graphs among the randomly generated graphs
    connected = 0;
    
    parpool(24);
    parfor iter=1:n_iter
        try 
            n = 5 + randi(15);
    
            % Construct adjacency matrices
            A1 = randi(3, n, n) - 1;
            A1 = triu(A1,1) + diag(diag(A1)) + triu(A1,1).';
            
            % Random permutation
            perm = randperm(n);
            A2 = A1(perm,perm);

            % Get transition matrices
            P1 = A1 ./ sum(A1, 2);
            P2 = A2 ./ sum(A2, 2);

            % Get distributions
            stat_dist1 = approx_stat_dist(P1, 100)';
            stat_dist2 = approx_stat_dist(P2, 100)';
            unif_dist1 = ones(n,1)/n;
            unif_dist2 = ones(n,1)/n;

            % Get cost function
            c = get_degree_cost(A1, A2);
    
            % Run algorithm
            [~, ~, otc_alignment] = exact_otc(P1, P2, c);
            [otsd_alignment, ~] = computeot_lp(c', stat_dist1, stat_dist2');
            otsd_alignment = reshape(otsd_alignment, n, n)';
            [~, gw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 1);
            [~, fgw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 0.5);

            % Get alignment
            [~,idx_otc] = max(otc_alignment,[],2);
            [~,idx_otsd] = max(otsd_alignment,[],2);
            [~,idx_gw] = max(gw_alignment,[],2);
            [~,idx_fgw] = max(fgw_alignment,[],2);

            % Check success
            otc_success = otc_success + check_isomorphism(idx_otc, A1, A2);
            otsd_success = otsd_success + check_isomorphism(idx_otsd, A1, A2);
            gw_success = gw_success + check_isomorphism(idx_gw, A1, A2);
            fgw_success = fgw_success + check_isomorphism(idx_fgw, A1, A2);
            
            % Update connected graph number
            connected = connected + 1;
        catch
        end
    end
    delete(gcp('nocreate'));
    
    % Display accuracy
    disp(['OTC mean accuracy: ', num2str(otc_success/connected)]);
    disp(['OTSD mean accuracy: ', num2str(otsd_success/connected)]);
    disp(['GW mean accuracy: ', num2str(gw_success/connected)]);
    disp(['FGW mean accuracy: ', num2str(fgw_success/connected)]);
    
    disp(connected);
    
    
    % 9.
    % Lollipop graph
    disp('Lollipop graph start')

    % Number of success for each algorithm
    otc_success = 0;
    otsd_success = 0;
    gw_success = 0;
    fgw_success = 0;
    
    % Number of connected graphs among the randomly generated graphs
    connected = 0;
    
    parpool(24);
    parfor iter=1:n_iter
        try 
            n1 = 6+randi(9); % number of nodes in the candy part
            n2 = 6+randi(9); % number of nodes in the stick part
            n = n1+n2;
            p = 0.5;

            % Construct adjacency matrices
            [~,A1] = stochastic_block_model([n1], [p]);
            A1 = triu(A1,1);
            for row=1:(n1-1)
                if A1(row,row+1) == 0
                    A1(row,row+1) = 1;
                end
            end
            A1(1,n1) = 1;
            A1 = A1 + A1.';
            
            A2 = zeros(n2,n2);
            for row=1:(n2-1)
                A2(row,row+1) = 1;
            end
            A2 = A2 + A2.';
            
            A = zeros(n,n);
            A(1:n1,1:n1) = A1;
            A(n1+1:n,n1+1:n) = A2;
            A(n1,n1+1) = 1;
            A(n1+1,n1) = 1;
            A1 = A;
            
            % Random permutation
            perm = randperm(n);
            A2 = A1(perm,perm);

            % Get transition matrices
            P1 = A1 ./ sum(A1, 2);
            P2 = A2 ./ sum(A2, 2);

            % Get distributions
            stat_dist1 = approx_stat_dist(P1, 100)';
            stat_dist2 = approx_stat_dist(P2, 100)';
            unif_dist1 = ones(n,1)/n;
            unif_dist2 = ones(n,1)/n;

            % Get cost function
            c = get_degree_cost(A1, A2);
    
            % Run algorithm
            [~, ~, otc_alignment] = exact_otc(P1, P2, c);
            [otsd_alignment, ~] = computeot_lp(c', stat_dist1, stat_dist2');
            otsd_alignment = reshape(otsd_alignment, n, n)';
            [~, gw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 1);
            [~, fgw_alignment] = fgw_dist(c, A1, A2, unif_dist1, unif_dist2, 1, 0.5);

            % Get alignment
            [~,idx_otc] = max(otc_alignment,[],2);
            [~,idx_otsd] = max(otsd_alignment,[],2);
            [~,idx_gw] = max(gw_alignment,[],2);
            [~,idx_fgw] = max(fgw_alignment,[],2);

            % Check success
            otc_success = otc_success + check_isomorphism(idx_otc, A1, A2);
            otsd_success = otsd_success + check_isomorphism(idx_otsd, A1, A2);
            gw_success = gw_success + check_isomorphism(idx_gw, A1, A2);
            fgw_success = fgw_success + check_isomorphism(idx_fgw, A1, A2);
            
            % Update connected graph number
            connected = connected + 1;
        catch
        end
    end
    delete(gcp('nocreate'));
    
    % Display accuracy
    disp(['OTC mean accuracy: ', num2str(otc_success/connected)]);
    disp(['OTSD mean accuracy: ', num2str(otsd_success/connected)]);
    disp(['GW mean accuracy: ', num2str(gw_success/connected)]);
    disp(['FGW mean accuracy: ', num2str(fgw_success/connected)]);

    disp(connected);

end