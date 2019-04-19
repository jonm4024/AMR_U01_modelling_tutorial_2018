% This script uses the virulence phenotypes measured for 
% 31 strains of Clostridium difficile from the paper
% Lewis et al 2017 (DOI: 10.1128/mBio.00885-17)

% data included in subdirectory 'data' %%%%%%%%%%%%%%%%%%%%%%%%%
%
% c.difficile.clinical.phenotypes.with.explanations.170606.xlsx:
% the phenotypes measured in Lewis et al 2017
%
% all32Strains_Patric.txt: 
% genes in iCN900 model present in each strain						
%
% growthCap.csv: 
% growth capabilities predicted by model for each strain
%
% cdiff_630_reference_NCBI.fasta:  
% genome of strain 630 used in iCN900

% reproduce the color map of C. diff the clades in Lewis et al
cMapClade = table({'Clade 1' 'Clade 2' 'NA' 'Clade 4' 'Clade 5'}',[0.19 0.25 0.52; 0.7 0.2 0.15;...
    0 0 0; 0.47 0.01 0.52; 0.22 0.39 0.12]);
cMapClade.Properties.VariableNames = {'clade', 'color'};

figure(1)
imagesc((1:height(cMapClade))')
set(gca, 'YTick', 1:height(cMapClade), 'YTickLabel', cMapClade.clade);
set(gca, 'XTick', []);
colormap(cMapClade.color)
title('C.diff clade color from Lewis et al 2017')

%% Load the virulence phenotypes from Lewis et al 2017 and cleanup data
% Read the virulence phenotypes
tVirulencePhenotypes = readtable('data/c.difficile.clinical.phenotypes.with.explanations.170606.xlsx',...
    'Range','A17:Q51' );
% replace percent of mice survived by percent of mice died;
tVirulencePhenotypes.mortality = 100 - tVirulencePhenotypes.survival;
tVirulencePhenotypes.survival = [];
% correct minor discrepencies with strain naming
tVirulencePhenotypes.Strain = upper(tVirulencePhenotypes.Strain);
tVirulencePhenotypes.Strain = strrep(tVirulencePhenotypes.Strain, 'WUP', 'WU');
% add the strain names as row indices
tVirulencePhenotypes.Properties.RowNames = tVirulencePhenotypes.Strain;

head(tVirulencePhenotypes)

%% read the tables with presence/absence of metabolic genes
tMetabolicGenesAll = readtable('data/all32Strains_Patric.txt');
% get name of strains to add as column names
strainsWithGenome = tMetabolicGenesAll.Properties.VariableNames(2:end);
% cleanup
strainsWithGenome = upper(strainsWithGenome);
strainsWithGenome = strrep(strainsWithGenome, 'WUP', 'WU');
% add strain names as column names
tMetabolicGenesAll.Properties.VariableNames(2:end) = strainsWithGenome; 

head(tMetabolicGenesAll)


%% Get predictions of growth capabilities obtained from model
tGrowthCapacities = readtable('data/growthCap.csv');
tGrowthCapacities.Properties.VariableNames{1} = 'Strain';
tGrowthCapacities.Strain = upper(tGrowthCapacities.Strain);
tGrowthCapacities.Strain = strrep(tGrowthCapacities.Strain, 'WUP', 'WU');

% select only the carbon sources where there is interesting variability
tGrowthSims = tGrowthCapacities(:, {'Strain' 'EX_acgam_e' 'EX_tre_e' 'EX_lac__L_e'});
tGrowthSims.Properties.RowNames = tGrowthSims.Strain;
tGrowthSims = tGrowthSims(strainsWithGenome, :);

% Add the predicted growth phenotypes to the virulence phenotypes matrix
tVirulencePhenotypes = tVirulencePhenotypes(ismember(tVirulencePhenotypes.Strain, strainsWithGenome), :);
tVirulencePhenotypes = tVirulencePhenotypes(strainsWithGenome, :);
tVirulencePhenotypes = [tVirulencePhenotypes tGrowthSims(:, 2:end)];


%% unsupoervized analyzis of virulence phenotypes and growth capacities
% predicted by model
phenotypesToAnalyze = {'patientSeverity' 'toxin_d7', 'spor_10000', 'vir_d2', 'vir_d3',...
    'vir_d6', 'vir_d14', 'lca_7hr', 'dca_7hr', 'mortality', 'EX_acgam_e', ...
    'EX_tre_e', 'EX_lac__L_e'};

% make a clustergram of the phenotypes
cgo = clustergram(tVirulencePhenotypes{:, phenotypesToAnalyze},...
    'RowLabels', tVirulencePhenotypes.Strain,...
    'Standardize', 'column',...
    'ColumnLabels', phenotypesToAnalyze,...
    'ColumnPDist', 'cosine', 'RowPDist', 'cosine')
close hidden;

% plot the clustergram
figure(3)
h = plot(cgo, gcf); 
set(h,'TickLabelInterpreter','none');
s = h.YTickLabel;
hCgo = gca;

% plot the color code of clades
figure(4)
set(gcf, 'Position', [1559, 936, 261, 317]);
s = s(end:-1:1);
g = [tVirulencePhenotypes.Clade(s)];
imagesc(g)
set(gca, 'YTick', 1:length(s), 'YTickLabel', s)
set(gca, 'XTick', 1:length(s), 'XTickLabel', {'Clade' 'MLST'})
colormap(cMapClade.color)

%% Presence/absence of metabolic genes shows no difference between X186A and WU4 
% It does show differences between WU8, a clade 2 isolate that caused severe 
% infection in the patient but not in mouse.
clade2Strains = tVirulencePhenotypes.Strain(tVirulencePhenotypes.Clade == 2);

% find all genes different betweeen the strains
m = sum(tMetabolicGenesAll{:, clade2Strains} > 0.95, 2);
idxDiffs = find(m ~= 0 & m~= length(clade2Strains));

tClade2 = tMetabolicGenesAll(idxDiffs, [{'locus'}; clade2Strains]);

% find the anotations and add it to table tClade2 for comparison
[genomeHeaders, ~] = fastaread('data/cdiff_630_reference_NCBI.fasta');

for i = 1:height(tClade2)
    idx = find(contains(genomeHeaders, [tClade2.locus{i} ' ']));
    tClade2.header{i} = genomeHeaders{idx};
end

% list of genes present in every Clade 2 strain except X186A and WU8
disp('Genes present in every Clade 2 strain except X186A and WU8:')
m = double(tClade2{:, 2:7}>0.95);
tClade2(~any(m(:, [1 4]), 2) & all(m(:, [2 3 5 6]), 2), :)

% list of genes present in every Clade 2 strain except X186A 
disp('Genes present in every Clade 2 strain except X186A:')
m = double(tClade2{:, 2:7}>0.95);
tClade2(~any(m(:, [4]), 2) & all(m(:, [1 2 3 5 6]), 2), :)


%% unsupoervized analyzis of virulence phenotypes and growth capacities (PCA)
[coeff,virulencePC,latent,tsquared,explained] = pca(normalize(tVirulencePhenotypes{:, phenotypesToAnalyze}),...
   'Centered',true);

figure(5)
subplot(1, 2, 1)
gscatter(virulencePC(:, 1), virulencePC(:, 2), tVirulencePhenotypes.Clade, cMapClade.color([1 2 4 5], :), [], 40)
colormap jet;
axis square
xlabel(sprintf('PC1 (%0.1f %% EV)', explained(1)))
ylabel(sprintf('PC2 (%0.1f %% EV)', explained(2)))
h = legend;
set(h, 'Location', 'East');
title('PCA of C. diff traits')

% biplot to explain the PCA
subplot(1, 2, 2)
h = biplot(coeff(:, 1:2),'Scores',virulencePC(:,1:2),'Varlabels',phenotypesToAnalyze);
axis square
set(gca,'TickLabelInterpreter','none');
set(gcf,'defaultTextInterpreter','none');
for i = 1:length(h)
    try
        set(h(i), 'Interpreter', 'none')
    catch
    end
end
title('biplot explaining the PCs')
xlabel(sprintf('PC1 (%0.1f %% EV)', explained(1)))
ylabel(sprintf('PC2 (%0.1f %% EV)', explained(2)))


%% Supervised analysis - search for phenotypes that explain virulence

% exclude only mouse mortality from the set of features to analyze
iPhenotypes = setdiff(phenotypesToAnalyze, {'mortality' 'patientSeverity'});

% use quadratic programming feature selection
classes = tVirulencePhenotypes.patientSeverity;
features = tVirulencePhenotypes{:, iPhenotypes};
featureScore = qpFeatureSelection(classes, features);

qpFeatures = table(iPhenotypes', featureScore);

% fit a predictor of virulence based on toxid_d7 alone
mdl1 = fitglm(tVirulencePhenotypes, 'patientSeverity ~ toxin_d7',...
    'Distribution','binomial','link','logit')

% fit a predictor of virulence based on the best features identified
mdl2 = fitglm(tVirulencePhenotypes, 'PredictorVars',qpFeatures{qpFeatures.featureScore>0.1, 'Var1'}, ...
    'ResponseVar','patientSeverity', ...
    'Distribution','binomial','link','logit')

% plot ROC and determine area under the cuve to 
figure(6)
set(gcf, 'Position', [412   222   564   583]);

[X,virulencePC,T,AUC1] = perfcurve(classes,mdl1.Fitted.Probability, 1);
plot(X,virulencePC, 'g-', 'LineWidth', 4)
hold on
[X,virulencePC,T,AUC2] = perfcurve(classes,mdl2.Fitted.Probability, 1);
plot(X,virulencePC, 'r-', 'LineWidth', 4)
plot([0 1],[0 1], 'k--')
hold off
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC for Classification by Logistic Regression')
axis equal tight
legend({char(mdl1.Formula) char(mdl2.Formula)},...
    'Location', 'southoutside', 'Interpreter', 'none')

fprintf('AUC for toxix_d7: %f\n', AUC1)
fprintf('AUC for model with best features: %f\n', AUC2)

%% try also which LASSO agrees with the result obtained by qpFeatureSelection
[B,FitInfo] = lassoglm(features,classes, 'binomial',  'CV',5, 'MCReps',10);
lassoPlot(B,FitInfo,'PlotType','CV');
qpFeatures.lassoFit = B(:, FitInfo.Index1SE)

%% find the genes unique to the strains that grow on trehalose
iTre = tVirulencePhenotypes.EX_tre_e > 0;
strainsTre = tVirulencePhenotypes.Strain(iTre);
strainsNotTre = tVirulencePhenotypes.Strain(~iTre);
iTreGenes = tMetabolicGenesAll.locus(all(tMetabolicGenesAll{:, strainsTre} > 0, 2) &...
    all(tMetabolicGenesAll{:, strainsNotTre} == 0, 2));

disp('Genes unique to strains predicted to grow on trehalose:')
for i = 1:length(iTreGenes)
    idx = find(contains(genomeHeaders, [iTreGenes{i} ' ']));
    disp(genomeHeaders{idx});
end

