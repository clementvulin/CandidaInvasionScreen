%% Loading data
load('/Users/cvulin/Documents/Zinkernagel/PCA_analysis/DataUsedForPCA_downstream.mat','DataPCA');

%%

rng(42); % For reproducibility
numBootstraps = 1000; % Number of bootstrap iterations
numClusters = 3; % Adjust based on data
numComponents = 2; % PCA components
nSamples = size(DataPCA, 1);

% Perform PCA on original dataset
[coeff, score] = pca(DataPCA);
reducedData = score(:, 1:numComponents);

% Cluster original PCA data
[idxOriginal, C] = kmeans(reducedData, numClusters, 'Replicates', 10);

% Calculate initial and witin group distances
    % Compute mean distances (centroid-based)
    withinDistOldOri = computeWithinGroupDistMean(reducedData, idxOriginal, numClusters);
    betweenDistOldOri = computeBetweenGroupDistMean(C);

    % Compute all distances (pairwise-based)
    withinDistNewOri = computeWithinGroupDistPairs(reducedData, idxOriginal, numClusters);
    betweenDistNewOri = computeBetweenGroupDistPairs(reducedData, idxOriginal, numClusters);


% Initialize storage for distances
withinDistOld = zeros(numBootstraps, 1);
betweenDistOld = zeros(numBootstraps, 1);
withinDistNew = zeros(numBootstraps, 1);
betweenDistNew = zeros(numBootstraps, 1);

for b = 1:numBootstraps
    % Bootstrap sampling
    sampleIdx = randsample(nSamples, nSamples, true); % Resample with replacement
    Data_bootstrap = DataPCA(sampleIdx, :);

    % Perform PCA on bootstrap sample
    [coeff_boot, score_boot] = pca(Data_bootstrap);
    reducedBootData = score_boot(:, 1:numComponents);

    % Cluster bootstrapped PCA data
    [idxBoot, C_boot] = kmeans(reducedBootData, numClusters, 'Replicates', 10);

%     % Compute pairwise distances in bootstrap sample
%     X_bootstrap = score_boot(idx, :); % Get resampled PCA scores
%     bootstrap_dist = pdist2(X_bootstrap, X_bootstrap);
% 
%     % Map bootstrap distances back to original sample positions
%     bootstrap_dists(idx, idx, b) = bootstrap_dist;


    % Compute mean distances (centroid-based)
    withinDistOld(b) = computeWithinGroupDistMean(reducedBootData, idxBoot, numClusters); %numCluster was C_boot
    betweenDistOld(b) = computeBetweenGroupDistMean(C_boot);

    % Compute all distances (pairwise-based)
    withinDistNew(b) = computeWithinGroupDistPairs(reducedBootData, idxBoot, numClusters);
    betweenDistNew(b) = computeBetweenGroupDistPairs(reducedBootData, idxBoot, numClusters);
end

%% Visualization of values
figure;
subplot(1,2,1);
boxplot([withinDistOld, withinDistNew], 'Labels', {'Old Method', 'New Method'});
title('Within-Group Distance Comparison');
ylabel('Distance');

subplot(1,2,2);
boxplot([betweenDistOld, betweenDistNew], 'Labels', {'Old Method', 'New Method'});
title('Between-Group Distance Comparison');
ylabel('Distance');


%% Visualization of ratios

% Compute the ratio of within-group to between-group distances
distanceRatioOld = withinDistOld ./ betweenDistOld;
distanceRatioNew = withinDistNew ./ betweenDistNew;

distanceRatioOldOri = withinDistOldOri/betweenDistOldOri;
distanceRatioNewOri = withinDistNewOri/betweenDistNewOri;
    
figure; tiledlayout(1,2); nexttile;hold on;
histogram(distanceRatioOld, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'FaceColor', 'b');
xline(distanceRatioOldOri, 'b', 'LineWidth', 2); % Mark original ratio
legend('Centroid-Based, bootstrap', 'Centroid-Based, original');
xlabel('Within-Group / Between-Group Distance Ratio');
ylabel('Density');
title('Within/between distance ratios, centroid');

nexttile; hold on;
histogram(distanceRatioNew, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'FaceColor', 'r');
xline(distanceRatioNewOri, 'r', 'LineWidth', 2); % Mark original ratio
legend('Pair-Based, bootstrap', 'Pair-Based, original');
xlabel('Within-Group / Between-Group Distance Ratio');
ylabel('Density');
title('Within/between distance ratios, pair');


grid on;
hold off;

%%
computeBetweenGroupDistMean(C)

%% functions

% Function to compute within-group distance
function withinDist = computeWithinGroupDistMean(data, labels, numClusters)
    withinDist = 0;
    for i = 1:numClusters
        clusterPoints = data(labels == i, :);
        clusterMean = mean(clusterPoints, 1);
        withinDist = withinDist + sum(vecnorm(clusterPoints - clusterMean, 2, 2).^2);
    end
    withinDist = withinDist / size(data, 1); % Normalize by number of points
end

% Function to compute between-group distance
function betweenDist = computeBetweenGroupDistMean(clusterCenters)
    betweenDist = 0;
    numClusters = size(clusterCenters, 1);
    for i = 1:numClusters
        for j = i+1:numClusters
            betweenDist = betweenDist + norm(clusterCenters(i, :) - clusterCenters(j, :))^2;
        end
    end
    betweenDist = betweenDist / (numClusters * (numClusters - 1) / 2); % Normalize by pairs
end


% Function to compute within-group distance using pairwise distances
function withinDist = computeWithinGroupDistPairs(data, labels, numClusters)
    withinDist = 0;
    totalPairs = 0;
    
    for i = 1:numClusters
        clusterPoints = data(labels == i, :);
        numPoints = size(clusterPoints, 1);
        
        if numPoints > 1
            % Compute pairwise distances
            D = pdist(clusterPoints, 'euclidean');
            withinDist = withinDist + sum(D.^2); % Sum of squared distances
            totalPairs = totalPairs + length(D); % Count pairs
        end
    end
    
    % Normalize by the number of point pairs
    if totalPairs > 0
        withinDist = withinDist / totalPairs;
    end
end

% Function to compute between-group distance using pairwise distances
function betweenDist = computeBetweenGroupDistPairs(data, labels, numClusters)
    betweenDist = 0;
    totalPairs = 0;
    
    for i = 1:numClusters
        for j = i+1:numClusters
            points_i = data(labels == i, :);
            points_j = data(labels == j, :);
            
            if ~isempty(points_i) && ~isempty(points_j)
                % Compute pairwise distances between all points in different clusters
                D = pdist2(points_i, points_j, 'euclidean');
                betweenDist = betweenDist + sum(D(:).^2); % Sum of squared distances
                totalPairs = totalPairs + numel(D); % Count pairs
            end
        end
    end
    
    % Normalize by the number of point pairs
    if totalPairs > 0
        betweenDist = betweenDist / totalPairs;
    end
end
