%%
% classification_experiment.m
%
% Nearest neighbor classification experiment via graph optimal transport.
%%

function classification_experiment(dataset, longleaf, testing)
    disp(dataset);

    % Setup directories
    rng(2021);
    expid = strrep(num2str(now, '%f'), '.', '_');
    format longG;
    if longleaf
        addpath('/nas/longleaf/home/koconn/GraphOTC/GraphOTC/', '/nas/longleaf/home/koconn/GraphOTC/GraphOTC/otc/', '/nas/longleaf/home/koconn/GraphOTC/GraphOTC/otc/ot_algorithms/', '/nas/longleaf/home/koconn/GraphOTC/GraphOTC/utils/', '/nas/longleaf/home/koconn/GraphOTC/GraphOTC/experiments/');
        datadir = ['/pine/scr/k/o/koconn/GraphOTC/Data/' dataset '/'];
        expdir = [datadir expid '/'];
    else
        datadir = ['C:\Users\oconn\Documents\Research\OptimalJoinings\GraphOTC\Data\' dataset '\'];
        expdir = [datadir expid '\'];
    end
    mkdir(expdir);

    % Experiment parameters
    if testing
        runonsubset = 1;
        subset_size = 30;
        n_iter = 1;
    else
        runonsubset = 1;
        subset_size = 250;
        n_iter = 5;
    end
        
    % Algorithm switches
    otc = 1;
    fgw = 1;
    gw = 1;
    copt = 1;
    otsd = 1;

    % Load data
    data_adj = readtable([datadir dataset '_A.txt'], 'ReadVariableNames', 0);
    data_node_label = readtable([datadir dataset '_node_labels.txt'], 'ReadVariableNames', 0);
    data_graph_indicator = readtable([datadir dataset '_graph_indicator.txt'], 'ReadVariableNames', 0);
    data_graph_labels = readtable([datadir dataset '_graph_labels.txt'], 'ReadVariableNames', 0);
    data_node_label = table2array(data_node_label);

    % Split into test and train
    n_graphs = size(data_graph_labels, 1);
    [train_idxs_mat, test_idxs_mat] = get_n_training_sets(n_graphs, n_iter);
    if runonsubset && (n_graphs > subset_size)
        for i=1:n_iter
            perm = randperm(n_graphs);
            test_idxs_mat(:,i) = perm(test_idxs_mat(:,i));
            train_idxs_mat(:,i) = perm(train_idxs_mat(:,i));
        end
        test_idxs_mat = test_idxs_mat(1:(subset_size-round(0.8*subset_size)),:);
        train_idxs_mat = train_idxs_mat(1:round(0.8*subset_size),:);
    end
    n_test = size(test_idxs_mat, 1);
    n_train = size(train_idxs_mat, 1);

    otc_dist_tensor = 1e5*ones(n_iter, n_train, n_test);
    gw_dist_tensor = 1e5*ones(n_iter, n_train, n_test);
    fgw_dist_tensor = 1e5*ones(n_iter, n_train, n_test);
    otsd_dist_tensor = 1e5*ones(n_iter, n_train, n_test);
    copt_dist_tensor = 1e5*ones(n_iter, n_train, n_test);
    
    data_adj = unique(sort(table2array(data_adj), 2), 'rows');
    full_graph = graph(data_adj(:,1), data_adj(:,2));

    %parpool('local', feature('numcores'));
    %parfor i=1:n_train    
    for iter=1:n_iter
        disp(['Iteration ' num2str(iter)]);

        train_idxs = train_idxs_mat(:,iter);
        test_idxs = test_idxs_mat(:,iter);
        train_labels = table2array(data_graph_labels(train_idxs,1));
        test_labels = table2array(data_graph_labels(test_idxs,1));

        % Make distance matrices
        otc_dist_mat = 1e5*ones(n_train, n_test);
        gw_dist_mat = 1e5*ones(n_train, n_test);
        fgw_dist_mat = 1e5*ones(n_train, n_test);
        otsd_dist_mat = 1e5*ones(n_train, n_test);
        copt_dist_mat = 1e5*ones(n_train, n_test);
            
        % Run on test set
        parpool('local', feature('numcores'));
        parfor i=1:n_train
            disp(['Training sample ' num2str(i)]);

            % Construct graph i from training set
            g1_nodes = find(table2array(data_graph_indicator) == train_idxs(i));
            g1 = subgraph(full_graph, g1_nodes);
            A1 = full(adjacency(g1));
            n1 = size(A1, 1);

            % Construct transition matrix making sure it is ergodic
            P1 = adj_to_trans(A1+1e-2);
            stat_dist1 = get_best_stat_dist(P1, normrnd(0, 1, n1, 1));

            for j=1:n_test
                % Construct graphs
                g2_nodes = find(table2array(data_graph_indicator) == test_idxs(j));
                g2 = subgraph(full_graph, g2_nodes);
                A2 = full(adjacency(g2));
                n2 = size(A2, 1);

                % Construct transition matrix making sure it is ergodic
                P2 = adj_to_trans(A2+1e-2);
                stat_dist2 = get_best_stat_dist(P2, normrnd(0, 1, n2, 1));  

                % Compute cost between nodes
                c = zeros(n1, n2);
                if dataset == "Cuneiform"
                    for l=1:n1
                        for k=1:n2
                            c(l,k) = (data_node_label(g1_nodes(l),1) ~= data_node_label(g2_nodes(k),1)) | (data_node_label(g1_nodes(l),2) ~= data_node_label(g2_nodes(k),2));
                        end
                    end
                else
                    for l=1:n1
                        for k=1:n2
                            c(l,k) = (data_node_label(g1_nodes(l),:) ~= data_node_label(g2_nodes(k),:));
                        end
                    end
                end
                    
                % Compute OT costs
                if otc
                    try 
                        otc_dist_mat(i,j) = entropic_otc(P1, P2, c, 25, 50, 100, 50, 0);
                    end
                end
                if gw
                    try
                        gw_dist_mat(i,j) = fgw_dist(c, A1, A2, stat_dist1, stat_dist2, 1, 1.0);
                    end
                end
                if fgw
                    try
                        fgw_dist_mat(i,j) = fgw_dist(c, A1, A2, stat_dist1, stat_dist2, 1, 0.5);
                    end
                end
                if otsd
                    try
                        [~, otsd_dist_mat(i,j)] = computeot_lp(c', stat_dist1, stat_dist2');
                    end
                end
                if copt
                    try
                        copt_dist_mat(i,j) = copt_dist(A1, A2, 10);
                    end
                end
            end
        end
        delete(gcp('nocreate'));

        otc_dist_tensor(iter,:,:) = otc_dist_mat;
        gw_dist_tensor(iter,:,:) = gw_dist_mat;
        fgw_dist_tensor(iter,:,:) = fgw_dist_mat;
        otsd_dist_tensor(iter,:,:) = otsd_dist_mat;
        copt_dist_tensor(iter,:,:) = copt_dist_mat;

        %%
        % Classify graphs
        k_vec = [5];
        otc_err = zeros(length(k_vec), 1);
        gw_err = zeros(length(k_vec), 1);
        fgw_err = zeros(length(k_vec), 1);
        otsd_err = zeros(length(k_vec), 1);
        copt_err = zeros(length(k_vec), 1);

        for j=1:length(k_vec)
            n_neigh = k_vec(j);

            % Perform classification
            otc_labels = knn(otc_dist_mat, train_labels, n_neigh);
            gw_labels = knn(gw_dist_mat, train_labels, n_neigh);
            fgw_labels = knn(fgw_dist_mat, train_labels, n_neigh);
            otsd_labels = knn(otsd_dist_mat, train_labels, n_neigh);
            copt_labels = knn(copt_dist_mat, train_labels, n_neigh);

            % Compute errors
            otc_err(j) = mean(otc_labels ~= test_labels);
            gw_err(j) = mean(gw_labels ~= test_labels);
            fgw_err(j) = mean(fgw_labels ~= test_labels);
            otsd_err(j) = mean(otsd_labels ~= test_labels);
            copt_err(j) = mean(copt_labels ~= test_labels);
        end

        % Write errors
        disp('Classification error');
        algorithms = {'OTC'; 'GW'; 'FGW'; 'OT_StatDist'; 'COPT'};
        algorithms_used = algorithms(find([otc gw fgw otsd copt]),1);
        errors = [otc_err gw_err fgw_err otsd_err copt_err];
        errors_computed = errors(:,find([otc gw fgw otsd copt]));
        errors_table = array2table(errors_computed', 'RowNames', algorithms_used);
        disp(errors_table);
        writetable(errors_table, [expdir dataset '_errors' num2str(iter) '.txt'], 'WriteRowNames', 1);

        % Write accuracies
        disp('Classification accuracy');
        accuracies = [1-otc_err 1-gw_err 1-fgw_err 1-otsd_err 1-copt_err];
        accuracies_computed = accuracies(:,find([otc gw fgw otsd copt]));
        accuracies_table = array2table(accuracies_computed', 'RowNames', algorithms_used);
        disp(accuracies_table);
        writetable(accuracies_table, [expdir dataset '_accuracies' num2str(iter) '.txt'], 'WriteRowNames', 1);
    end

    %% 
    % Store distance tensors
    save([expdir 'all_variables.mat']);
end