% PCA analysis

load('/Users/cvulin/Documents/Zinkernagel/PCA_analysis/DataUsedForPCA_downstream.mat','DataPCA');

% this was the initial anaylsis:

% basicPCAPlot(DataPCA)


% same with standardised data

% Standardize the data => this doesn't make sense since all in the same
% unit
% X_std = zscore(DataPCA);
% basicPCAPlot(X_std)

% basicPCAPlot(DataPCA)

% Correlation Matrix


% Check variance across variables
var_X = var(DataPCA');

% Compute covariance and correlation matrices
%cov_matrix = cov(X_std');
corrMatrix = corrcoef(DataPCA);

% Convert correlation to a distance metric (1 - absolute correlation)
distMatrix = 1 - abs(corrMatrix);

% Perform hierarchical clustering
Z = linkage(squareform(distMatrix), 'average'); % Use 'average', 'ward', or 'complete' as needed

% Get the ordering of variables from the clustering
order = optimalleaforder(Z, distMatrix);

% Reorder the correlation matrix
reorderedCorr = corrMatrix(order, order);

%% Visualize the correlation matrix
figure; imagesc(reorderedCorr);
colorbar;
title('Correlation Matrix Heatmap');

%% Optionally, plot the dendrogram
figure;
dendrogram(Z, 'Reorder', order);
title('Correlation-based Hierarchical Clustering');


%% remove variables with a corr thresh above 0.5 ?

thresh=linspace(0.5,0.999,100);

for tis=1:numel(thresh)
corrThreshold = thresh(tis);  % Set correlation cutoff
distThreshold = 1 - corrThreshold;  % Convert correlation to distance

% Cut the dendrogram at the corresponding distance
clusters = cluster(Z, 'cutoff', distThreshold, 'criterion', 'distance');

% Find Number of Clusters
numClusters(tis) = numel(unique(clusters));
end

figure; plot (thresh,numClusters,'LineWidth',2)
xlabel('Thrshold proposed');
ylabel( 'number of clusters');
set(gca, 'TickDir','out', 'FontSize',20, 'Box','off','LineWidth',2)

%% remove variables with a corr thresh above which value?

thresh=linspace(0,0.999,1000);

for tis=1:numel(thresh)
corrThreshold = thresh(tis);  % Set correlation cutoff
distThreshold = 1 - corrThreshold;  % Convert correlation to distance

% Cut the dendrogram at the corresponding distance
clusters = cluster(Z, 'cutoff', distThreshold, 'criterion', 'distance');

% Find Number of Clusters
numClusters(tis) = numel(unique(clusters));
end

figure; plot (thresh,numClusters,'LineWidth',2)
xlabel('Merging threshold');
ylabel( 'Number of clusters');
set(gca, 'TickDir','out', 'FontSize',20, 'Box','off','LineWidth',2)

xline(0.6, 'LineWidth',2)

%% create correlation map for the 0.6 threshold:

% Compute clusters
corrThreshold = 0.6;  
distThreshold = 1 - corrThreshold;  
clusters = cluster(Z, 'cutoff', distThreshold, 'criterion', 'distance');
numClusters = numel(unique(clusters));

% Select Representative Variables
selectedVars = [];
for i = 1:numClusters
    clusterVars = find(clusters == i);  
    avgCorr = mean(abs(corrMatrix(clusterVars, clusterVars)), 2); 
    [~, bestVarIdx] = max(avgCorr); 
    selectedVars = [selectedVars; clusterVars(bestVarIdx)];
end
selectedVars = find(ismember(order, selectedVars)); % Adjust for reordering

% Define the custom colormap with white at 0.6
nBlue = 80;  % Number of steps from -1 to 0.6
nRed = 20;   % Number of steps from 0.6 to 1
% Blue to White (for -1 to 0.6)
blueToWhite = [linspace(0,1,nBlue)', linspace(0,1,nBlue)', ones(nBlue,1)];
% White to Red (for 0.6 to 1)
whiteToRed = [ones(nRed,1), linspace(1,0,nRed)', linspace(1,0,nRed)'];
% Combine both parts
cmap = [blueToWhite; whiteToRed];

% Plot Heatmap
figure;
imagesc(reorderedCorr, [-1 1]); % Set color scale from -1 to 1
colormap(cmap);
colorbar;
title('Correlation Matrix Heatmap');
hold on;
caxis([-1 1]); % Ensures correct color mapping
set(gca, 'TickDir','out', 'FontSize',20, 'Box','off','LineWidth',2)

% Add Stars Above Selected Variables
for i = 1:length(selectedVars)
    text(selectedVars(i), 0, '*', 'FontSize', 14, 'Color', 'k', 'HorizontalAlignment', 'center');
end

% Draw Cluster Separation Lines
hold on;
edges = find(diff(clusters(order)) ~= 0); % Find boundaries between clusters
for i = 1:length(edges)
    x = edges(i) + 0.5; % Line position
    y = size(reorderedCorr, 1) + 0.5; 
    plot([x, x], [0.5, y], 'k-', 'LineWidth', 1.5); % Vertical line
    plot([0.5, y], [x, x], 'k-', 'LineWidth', 1.5); % Horizontal line
end

% Adjust Axes
xticks(1:length(order));
yticks(1:length(order));
xticklabels(order);
yticklabels(order);
axis square;
hold off;



%% Select Representative Variables from Each Cluster


corrThreshold = 0.6;  % Set correlation cutoff
distThreshold = 1 - corrThreshold;  % Convert correlation to distance

% Cut the dendrogram at the corresponding distance
clusters = cluster(Z, 'cutoff', distThreshold, 'criterion', 'distance');

% Find Number of Clusters
numClusters = numel(unique(clusters));


selectedVars = [];  % Store selected variable indices
for i = 1:numClusters
    clusterVars = find(clusters == i);  % Get all variables in this cluster
    avgCorr = mean(abs(corrMatrix(clusterVars, clusterVars)), 2); % Compute average correlation
    [~, bestVarIdx] = max(avgCorr); % Pick the variable with the highest average correlation
    selectedVars = [selectedVars; clusterVars(bestVarIdx)];
end

% Extract the reduced dataset
DataPCA_reduced = DataPCA(:, selectedVars);

% Step 5: Perform PCA on Original and Reduced Data
[coeffFull, scoreFull, latentFull] = pca(DataPCA); % Original PCA
[coeffRed, scoreRed, latentRed] = pca(DataPCA_reduced); % PCA on reduced data

% Step 6: Compare Explained Variance
explainedFull = (latentFull) / sum(latentFull);
explainedRed = (latentRed) / sum(latentRed);

figure;
plot(explainedFull, '-o', 'DisplayName', 'Original PCA');
hold on;
plot(explainedRed, '-s', 'DisplayName', 'Cluster-Reduced PCA');
xlabel('Number of Principal Components');
ylabel('Cumulative Explained Variance');
legend;
title('Comparison of PCA Explained Variance');
grid on;

% 
basicPCAPlot(DataPCA)
basicPCAPlot(DataPCA_reduced)

%% How the PCAs look like

% Select Representative Variables from Each Cluster

for corrThreshold = 0.2:0.1:0.4  % Set correlation cutoff

distThreshold = 1 - corrThreshold;  % Convert correlation to distance

% Cut the dendrogram at the corresponding distance
clusters = cluster(Z, 'cutoff', distThreshold, 'criterion', 'distance');

% Find Number of Clusters
numClusters = numel(unique(clusters));


selectedVars = [];  % Store selected variable indices
for i = 1:numClusters
    clusterVars = find(clusters == i);  % Get all variables in this cluster
    avgCorr = mean(abs(corrMatrix(clusterVars, clusterVars)), 2); % Compute average correlation
    [~, bestVarIdx] = max(avgCorr); % Pick the variable with the highest average correlation
    selectedVars = [selectedVars; clusterVars(bestVarIdx)];
end

% Extract the reduced dataset
DataPCA_reduced = DataPCA(:, selectedVars);

basicPCAPlot(DataPCA_reduced);
title (['correlation Threshold: ', num2str(corrThreshold)])

end

% Extract the reduced dataset
DataPCA_reduced = DataPCA(:, selectedVars);

basicPCAPlot(DataPCA);
title ('original PCA')


%%
basicPCAPlot(DataPCA);

%% Compare PCA Sample Grouping Using Within/Between Distances

% Step 7.1: Define 3 groups based on original PCA
numGroups = 3;
groupLabels = kmeans(scoreFull(:, 1:2), numGroups); % Cluster samples using first 2 PCs

% Step 7.2: Compute cluster centers for both PCA spaces
clusterCentersFull = zeros(numGroups, 2);
clusterCentersRed = zeros(numGroups, 2);

for i = 1:numGroups
    clusterCentersFull(i, :) = mean(scoreFull(groupLabels == i, 1:2), 1);
    clusterCentersRed(i, :) = mean(scoreRed(groupLabels == i, 1:2), 1);
end

% Step 7.3: Compute within-group and between-group distances
within_full = computeWithinGroupDist(scoreFull(:, 1:2), groupLabels, numGroups);
between_full = computeBetweenGroupDist(clusterCentersFull);

within_red = computeWithinGroupDist(scoreRed(:, 1:2), groupLabels, numGroups);
between_red = computeBetweenGroupDist(clusterCentersRed);

% Step 7.4: Compute within/between ratio
ratio_full = within_full / between_full;
ratio_red = within_red / between_red;

% Step 7.5: Plot the ratios for comparison
figure; tiledlayout(2,1): nexttile
bar([ratio_full, ratio_red]);
set(gca, 'XTickLabel', {'Original PCA', 'Reduced PCA'});
ylabel('Within-Group / Between-Group Distance Ratio');
title('Group Distance Preservation, centroid');
grid on;

% Step 7.6: Display results
fprintf('Original PCA: Within/Between Distance Ratio = %.4f\n', ratio_full);
fprintf('Reduced PCA: Within/Between Distance Ratio = %.4f\n', ratio_red);


%%
%Step 7: Compare PCA Sample Grouping Using Pairwise Distances

% Step 7.1: Define 3 groups based on original PCA
numGroups = 3;
groupLabels = kmeans(scoreFull(:, 1:2), numGroups); % Cluster samples using first 2 PCs

% Step 7.2: Compute within-group and between-group distances using pairwise method
within_full = computeWithinGroupDistPairs(scoreFull(:, 1:2), groupLabels, numGroups);
between_full = computeBetweenGroupDistPairs(scoreFull(:, 1:2), groupLabels, numGroups);

within_red = computeWithinGroupDistPairs(scoreRed(:, 1:2), groupLabels, numGroups);
between_red = computeBetweenGroupDistPairs(scoreRed(:, 1:2), groupLabels, numGroups);

% Step 7.3: Compute within/between ratio
ratio_full = within_full / between_full;
ratio_red = within_red / between_red;

% Step 7.4: Plot the ratios for comparison
figure;
bar([ratio_full, ratio_red]);
set(gca, 'XTickLabel', {'Original PCA', 'Reduced PCA'});
ylabel('Within-Group / Between-Group Distance Ratio');
title('Comparison of PCA Group Structure Preservation (Pairwise Distances)');
grid on;

% Step 7.5: Display results
fprintf('Original PCA: Within/Between Distance Ratio = %.4f\n', ratio_full);
fprintf('Reduced PCA: Within/Between Distance Ratio = %.4f\n', ratio_red);


%% distance ratios using less variables
% Define Thresholds
corrThresholds = [0.5, 0.6, 0.7, 0.8, 0.9];
numThresholds = length(corrThresholds);

% Preallocate storage for results
ratios_centroid = zeros(numThresholds, 2); % [Original, Reduced]
ratios_pairwise = zeros(numThresholds, 2);

% Compute Original PCA and Fixed Clustering
[coeffFull, scoreFull, latentFull] = pca(DataPCA);
numGroups = 3;
groupLabels = kmeans(scoreFull(:, 1:2), numGroups, 'Replicates', 10); % Use replicates for stability

% Compute cluster centers for original PCA
clusterCentersFull = zeros(numGroups, 2);
for i = 1:numGroups
    clusterCentersFull(i, :) = mean(scoreFull(groupLabels == i, 1:2), 1);
end

% Compute distances for original PCA once
within_full_centroid = computeWithinGroupDist(scoreFull(:, 1:2), groupLabels, numGroups);
between_full_centroid = computeBetweenGroupDist(clusterCentersFull);
ratio_full_centroid = within_full_centroid / between_full_centroid;

within_full_pairwise = computeWithinGroupDistPairs(scoreFull(:, 1:2), groupLabels, numGroups);
between_full_pairwise = computeBetweenGroupDistPairs(scoreFull(:, 1:2), groupLabels, numGroups);
ratio_full_pairwise = within_full_pairwise / between_full_pairwise;

% Step 3: Loop Over Correlation Thresholds
for t = 1:numThresholds
    corrThreshold = corrThresholds(t);
    distThreshold = 1 - corrThreshold;

    % Select Representative Variables for Current Threshold
    clusters = cluster(Z, 'cutoff', distThreshold, 'criterion', 'distance');
    numClusters = numel(unique(clusters));

    selectedVars = [];
    for i = 1:numClusters
        clusterVars = find(clusters == i);
        avgCorr = mean(abs(corrMatrix(clusterVars, clusterVars)), 2);
        [~, bestVarIdx] = max(avgCorr);
        selectedVars = [selectedVars; clusterVars(bestVarIdx)];
    end
    numVarsKept(t)=numel(selectedVars);

    DataPCA_reduced = DataPCA(:, selectedVars);
    [coeffRed, scoreRed, latentRed] = pca(DataPCA_reduced);

    % Compute cluster centers for reduced PCA
    clusterCentersRed = zeros(numGroups, 2);
    for i = 1:numGroups
        clusterCentersRed(i, :) = mean(scoreRed(groupLabels == i, 1:2), 1);
    end

    % Compute distances for reduced PCA
    within_red_centroid = computeWithinGroupDist(scoreRed(:, 1:2), groupLabels, numGroups);
    between_red_centroid = computeBetweenGroupDist(clusterCentersRed);
    ratios_centroid(t, :) = [ratio_full_centroid, within_red_centroid / between_red_centroid];

    within_red_pairwise = computeWithinGroupDistPairs(scoreRed(:, 1:2), groupLabels, numGroups);
    between_red_pairwise = computeBetweenGroupDistPairs(scoreRed(:, 1:2), groupLabels, numGroups);
    ratios_pairwise(t, :) = [ratio_full_pairwise, within_red_pairwise / between_red_pairwise];

    fprintf('Threshold %.1f - Reduced PCA (Centroid): %.4f\n', corrThreshold, ratios_centroid(t, 2));
    fprintf('Threshold %.1f - Reduced PCA (Pairwise): %.4f\n', corrThreshold, ratios_pairwise(t, 2));
end

% Create figure with two subplots
figure;
tiledlayout(2,1);

% Centroid-based plot
nexttile;
bar(ratios_centroid);
hold on;
for i = 1:length(corrThresholds)
    text(i, ratios_centroid(i) + 0.05, sprintf('%d vars', numVarsKept(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'k');
end
hold off;
set(gca, 'XTickLabel', arrayfun(@num2str, corrThresholds, 'UniformOutput', false));
ylabel('Within-Group / Between-Group Distance Ratio');
title('Group Distance Preservation (Centroid)');
grid on;

% Pairwise-based plot
nexttile;
bar(ratios_pairwise);
hold on;
for i = 1:length(corrThresholds)
    text(i, ratios_pairwise(i) + 0.05, sprintf('%d vars', numVarsKept(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'k');
end
hold off;
set(gca, 'XTickLabel', arrayfun(@num2str, corrThresholds, 'UniformOutput', false));
ylabel('Within-Group / Between-Group Distance Ratio');
title('Group Distance Preservation (Pairwise)');
grid on;

%% Picking random sets of variables

% Step 1: Define Parameters
corrThreshold = 0.3;  % Set correlation cutoff
distThreshold = 1 - corrThreshold;  % Convert correlation to distance
numRuns = 1000;  % Number of random selections

% Cut the dendrogram at the corresponding distance
clusters = cluster(Z, 'cutoff', distThreshold, 'criterion', 'distance');

% Find Number of Clusters
numClusters = numel(unique(clusters));
numSamples = size(DataPCA, 1);
numPCs = 2; % Number of principal components to keep

% Perform PCA on the original dataset
[coeffFull, scoreFull, ~] = pca(DataPCA);

% Select two anchor samples (randomly chosen)
anchorIdx = [5,10];
anchorFull = scoreFull(anchorIdx, 1:numPCs); % Their projection in original PCA

% Initialize matrix to store cumulative aligned scores
scoreSum = zeros(numSamples, numPCs);
validRuns = 0;  % Track number of valid PCA runs

% Step 2: Run PCA 1000 Times with Random Variable Selection and Alignment
for run = 1:numRuns
    selectedVars = [];  % Store selected variable indices
    
    for i = 1:numClusters
        clusterVars = find(clusters == i);  % Get all variables in this cluster
        if ~isempty(clusterVars)
            randIdx = randi(length(clusterVars));  % Pick a random variable
            selectedVars = [selectedVars; clusterVars(randIdx)];
        end
    end
    
    % Perform PCA on reduced dataset
    DataPCA_reduced = DataPCA(:, selectedVars);
    [coeffRed, scoreRed, ~] = pca(DataPCA_reduced);
    
    % Ensure anchor points exist in reduced PCA
    if size(scoreRed, 1) < max(anchorIdx) || size(scoreRed, 2) < numPCs
        continue; % Skip this run if anchor points or PCs are missing
    end
    
% Select the two anchor samples from the reduced PCA space
anchorRed = scoreRed(anchorIdx, 1:numPCs); % 2 × numPCs

% Apply Procrustes only on the two anchor points
[~, ~, transform] = procrustes(anchorFull, anchorRed);

% Apply the transformation to all PCA-reduced points
%scoreAligned = scoreRed(:, 1:numPCs) * transform.T;  % Apply only rotation
scoreAligned = scoreRed(:, 1:numPCs) * transform.T * transform.b + transform.c(1, :);
    
    % Accumulate aligned PCA scores
    scoreSum = scoreSum + scoreAligned;
    validRuns = validRuns + 1;
end

% Step 3: Compute and Plot the Average PCA Projection
if validRuns > 0
    scoreAvg = scoreSum / validRuns; % Compute average
    
    % Define colors for specific samples
    color = repmat([0 0 0], numSamples, 1); % Default black
    for str = [1, 5, 6], color(str, :) = [1 0 0]; end  % Red for specific samples
    for str = [7, 9], color(str, :) = [0.3 1 0]; end   % Green for specific samples

    % Compute small shift for text labels
    decalage = range(scoreAvg(:, 1:2)) * 0.01;
    
    % Create figure
    figure;
    scatter(scoreAvg(:, 1), scoreAvg(:, 2), 100, color, 'filled'); hold on;
    scatter(scoreAvg(anchorIdx, 1), scoreAvg(anchorIdx, 2), 200, 'o', 'MarkerEdgeColor', [0 0 1], 'LineWidth', 2);
    
    % Add labels to points
    for stri = 1:numSamples
        text(scoreAvg(stri, 1) + decalage(1), scoreAvg(stri, 2) + decalage(2), num2str(stri), 'FontSize', 15);
    end
    
    % Set axis properties
    set(gca, 'TickDir', 'out', 'FontSize', 20, 'Box', 'off', 'LineWidth', 2);
    
    % Labels and title
    xlabel('PC1 (Aligned & Averaged)');
    ylabel('PC2 (Aligned & Averaged)');
    title('Anchored Average PCA Projection Over 1000 Random Variable Sets');
    
    % Grid
    grid on;
else
    disp('No valid PCA runs were completed. Check data consistency.');
end

%% PCA invariance

% Picking random sets of variables

% Step 1: Define Parameters
figure; tiledlayout(3,3);
for corrThreshold=0.2:0.1:0.9

distThreshold = 1 - corrThreshold;  % Convert correlation to distance
numRuns = 1000;  % Number of random selections

% Cut the dendrogram at the corresponding distance
clusters = cluster(Z, 'cutoff', distThreshold, 'criterion', 'distance');

% Find Number of Clusters
numClusters = numel(unique(clusters));
numSamples = size(DataPCA, 1);
numPCs = 2; % Number of principal components to keep

% Perform PCA on the original dataset
[coeffFull, scoreFull, ~] = pca(DataPCA);

% Select two anchor samples (randomly chosen)
anchorIdx = [5,10];
anchorFull = scoreFull(anchorIdx, 1:numPCs); % Their projection in original PCA

% Initialize matrix to store cumulative aligned scores
scoreSum = zeros(numSamples, numPCs);
validRuns = 0;  % Track number of valid PCA runs

% Step 2: Run PCA 1000 Times with Random Variable Selection and Alignment
for run = 1:numRuns
    selectedVars = [];  % Store selected variable indices
    
    for i = 1:numClusters
        clusterVars = find(clusters == i);  % Get all variables in this cluster
        if ~isempty(clusterVars)
            randIdx = randi(length(clusterVars));  % Pick a random variable
            selectedVars = [selectedVars; clusterVars(randIdx)];
        end
    end
    
    % Perform PCA on reduced dataset
    DataPCA_reduced = DataPCA(:, selectedVars);
    [coeffRed, scoreRed, ~] = pca(DataPCA_reduced);
    
    % Ensure anchor points exist in reduced PCA
    if size(scoreRed, 1) < max(anchorIdx) || size(scoreRed, 2) < numPCs
        continue; % Skip this run if anchor points or PCs are missing
    end
    
% Select the two anchor samples from the reduced PCA space
anchorRed = scoreRed(anchorIdx, 1:numPCs); % 2 × numPCs

% Apply Procrustes only on the two anchor points
[~, ~, transform] = procrustes(anchorFull, anchorRed);

% Apply the transformation to all PCA-reduced points
%scoreAligned = scoreRed(:, 1:numPCs) * transform.T;  % Apply only rotation
scoreAligned = scoreRed(:, 1:numPCs) * transform.T * transform.b + transform.c(1, :);
    
    % Accumulate aligned PCA scores
    scoreSum = scoreSum + scoreAligned;
    validRuns = validRuns + 1;
end

% Step 3: Compute and Plot the Average PCA Projection
if validRuns > 0
    scoreAvg = scoreSum / validRuns; % Compute average
    
    % Define colors for specific samples
    color = repmat([225 35 87]/255, numSamples, 1); % Default black
    for str = [1, 5, 6], color(str, :) = [64 117 211]/255; end  % Red for specific samples
    for str = [7, 9], color(str, :) = [64 35 211]/255; end   % Green for specific samples

    % Compute small shift for text labels
    decalage = range(scoreAvg(:, 1:2)) * 0.04;
    
    % Create figure
    nexttile;
    scatter(scoreAvg(:, 1), scoreAvg(:, 2), 100, color, 'filled'); hold on;
    scatter(scoreAvg(anchorIdx, 1), scoreAvg(anchorIdx, 2), 200, 'o', 'MarkerEdgeColor', [0 0 1], 'LineWidth', 2);
    
    % Add labels to points
    legs={'ReW','Th','Na','In1','ReO','In2','Wd','Tr','Ab','Lab'};
    for stri = 1:numSamples
        text(scoreAvg(stri, 1) + decalage(1), scoreAvg(stri, 2) + decalage(2), legs{stri}, 'FontSize', 15);
    end
    
    % Set axis properties
    set(gca, 'TickDir', 'out', 'FontSize', 20, 'Box', 'off', 'LineWidth', 2,'XTickLabel',[],'YTickLabel',[]);
    
    % Labels and title
    xlabel('PC1');
    ylabel('PC2');
    
    title(['T=',num2str(corrThreshold),', ', num2str(numClusters),'var.']);
    
else
    disp('No valid PCA runs were completed. Check data consistency.');
end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function basicPCAPlot(DataPCA)

[coeff,score,latent] = pca(DataPCA);

% labels for cond number
varlabel={};
for cond=1:size(DataPCA,2)
    varlabel{end+1}=['Cond',num2str(cond)]; 
end


figure; tiledlayout(1,3);nexttile;
obslabel={};for str=1:10; obslabel{end+1}=['s', num2str(str)]; end 
bar(latent./(sum(latent))); xlabel('Axis number'); ylabel('% Var'); set(gca, 'LineWidth',2,'FontSize',15, 'Box', 'off','TickDir','out');
set(gca, 'TickDir','out', 'FontSize',20, 'Box','off','LineWidth',2)

nexttile
biplot(coeff(:,1:2),'scores',score(:,1:2),'obslabels',obslabel,'varlabel', varlabel )
set(gca, 'TickDir','out', 'FontSize',20, 'Box','off','LineWidth',2)

nexttile
color=repmat([0 0 0],10,1);
for str=[1,5,6], color(str,:)=[1 0 0]; end
for str=[7,9], color(str,:)=[0.3 1 0]; end
scatter(score(:,1),score(:,2), 100,color, 'filled'); xline(0); yline(0); title('PCA on strains')
% add text
decalage=range(score(:,1:2))*0.01;
for stri=1:10
    text(score(stri,1)+decalage(1),score(stri,2)+decalage(2),num2str(stri))
end
set(gca, 'TickDir','out', 'FontSize',20, 'Box','off','LineWidth',2)

set(gcf, 'Position', [101         570        1376         420]);
end


% Function to compute within-group distance
function withinDist = computeWithinGroupDist(data, labels, numClusters)
    withinDist = 0;
    for i = 1:numClusters
        clusterPoints = data(labels == i, :);
        clusterMean = mean(clusterPoints, 1);
        withinDist = withinDist + sum(vecnorm(clusterPoints - clusterMean, 2, 2).^2);
    end
    withinDist = withinDist / size(data, 1); % Normalize by number of points
end

% Function to compute between-group distance
function betweenDist = computeBetweenGroupDist(clusterCenters)
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