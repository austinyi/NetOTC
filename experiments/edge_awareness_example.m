%%
% edge_awareness_example.m
%
%%

%% Setup graphs
G1 = graph([1 2 3 4 5 6 7 8],[2 3 4 5 6 7 8 1]);
d1 = zeros(8,2);
for i=1:8
    d1(i,1) = cos(pi/8+pi/4*(i-1));
    d1(i,2) = sin(pi/8+pi/4*(i-1));
end

G2 = graph([1 2 3 4 5 6 7],[2 3 4 5 6 7 8]);
d2 = d1;

G3 = graph([1 2 3 4 5 6 7],[2 3 4 5 6 7 8]);
d3 = zeros(8,2);
for i=1:8
    d3(i,1) = cos(pi/2+pi/7*(i-1));
    d3(i,2) = sin(pi/2+pi/7*(i-1));
end


%% Do NetOTC
% Get Adjacency matrix
A1 = adjacency(G1);
A1 = full(A1);

A2 = adjacency(G2);
A2 = full(A2);

A3 = adjacency(G3);
A3 = full(A3);


% Numbers of nodes of G1, G2, G3
n = 8;

% Get transition matrices
P1 = A1 ./ sum(A1, 2);
P2 = A2 ./ sum(A2, 2);
P3 = A3 ./ sum(A3, 2);

% Get distributions 
stat_dist1 = approx_stat_dist(P1, 100)';
stat_dist2 = approx_stat_dist(P2, 100)';
stat_dist3 = approx_stat_dist(P3, 100)';

unif_dist1 = ones(n,1)/n;
unif_dist2 = ones(n,1)/n;
unif_dist3 = ones(n,1)/n;

% Get cost matrix
c1 = zeros([n, n]);
for i=1:n
    for j=1:n
        c1(i, j) = sum((d2(i,:)-d1(j,:)).^2);
    end
end

c2 = zeros([n, n]);
for i=1:n
    for j=1:n
        c2(i, j) = sum((d2(i,:)-d3(j,:)).^2);
    end
end


% Run algorithm
% G2 vs G1
[otc_distance1, gotc1, otc_alignment1] = exact_otc(P2, P1, c1);
otc_distance1

[~, otsd_distance1] = computeot_lp(c1', stat_dist2, stat_dist1');
otsd_distance1

[fgw_distance1, fgw_alignment1] = fgw_dist(c1, A2, A1, unif_dist2, unif_dist1, 1, 0.5);
fgw_distance1

% G2 vs G3
[otc_distance2, gotc2, otc_alignment2] = exact_otc(P2, P3, c2);
otc_distance2

[~, otsd_distance2] = computeot_lp(c2', stat_dist2, stat_dist3');
otsd_distance2

[fgw_distance2, fgw_alignment2] = fgw_dist(c2, A2, A3, unif_dist2, unif_dist3, 1, 0.5);
fgw_distance2

%% Plot graphs G1, G2, G3
format longG;
savedir = ['/Users/bongsooyi/Documents/NetOTC/experiments/plot/'];
mkdir(savedir);

plot(G1,'XData',d1(:,1),'YData',d1(:,2),'LineWidth',6,'NodeFontSize',16) %'NodeColor','k','EdgeColor','k')
xlim([-1.2 1.2])
ylim([-1.2 1.2])
grid on
axis square
saveas(gcf,[savedir 'circle_G1.png'])


plot(G2,'XData',d2(:,1),'YData',d2(:,2),'LineWidth',6,'NodeFontSize',16) %'NodeColor','k','EdgeColor','k')
xlim([-1.2 1.2])
ylim([-1.2 1.2])
grid on
axis square
saveas(gcf,[savedir 'circle_G2.png'])

plot(G3,'XData',d3(:,1),'YData',d3(:,2),'LineWidth',6,'NodeFontSize',16)
xlim([-1.2 1.2])
ylim([-1.2 1.2])
grid on
axis square
saveas(gcf,[savedir 'circle_G3.png'])
