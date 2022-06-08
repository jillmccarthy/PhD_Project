%% %%%%%%% cTI %%%%%%%%%

% Extracts cTI stage (pseudotime), subtypes, and feature contributions for each subject from cTI output folders
% Calculates pearson correlations for cTI stage vs all relevant clinical/neuropsychological tests

clear
clc
close all

%set path to folder for full analysis
%path = '/data/SD_FTD/jill/MCM_cTI_files/cti/full'; 

%set path to folder for sensitivity analysis
path = '/data/SD_FTD/jill/MCM_cTI_files/cti/sensitivity_analysis';

% get all folder names (may include separate folders for combined analysis and for each individual modality depending on cTI analysis was done)
cti_folders = dir([path filesep '*_*']);

% set number of brain regions
N_regions = 78;

%read in list of brain region names
load('/data/SD_FTD/jill/MCM_files/region_names_78.mat')
region_names = region_names_78;

% set factor names
factor_names = ["GM density" "T1/T2 ratio" "fALFF" "FA" "MD"];

% read in genfi clinical table, for correlations between cTI and clinical/neuropsych tests, EYO
genfi_clinical = readtable('/export02/data/GENFI/spreadsheets/genfitable_clinical_neuropsych_bl.csv');

% remove values involving age for GRN199 - value is wrong!
genfi_clinical.AgeAtVisit(427) = 0;
genfi_clinical.AgeAtVisit(427) = NaN;
genfi_clinical.EYO(427) = NaN;
genfi_clinical.DiseaseDuration{427} = 'NA';

%read in clinical/eyo staging / genetic grouping file (to compare to cTI staging/subtyping)
clin_stages_table = readtable('/data/SD_FTD/jill/MCM_files/grouping_staging_4subtypes.txt', 'ReadVariableNames',false);

%%%% for sensitivity analysis, read in correct files (with only subset of subjects with all imaging modalities)

if contains(path, 'sensitivity_analysis')
    genfi_clinical = readtable('/data/SD_FTD/jill/MCM_files/genfi_clinical_sensitivity.txt');
    
    % remove values involving age for GRN199 - value is wrong!
    genfi_clinical.AgeAtVisit(170) = 0;
    genfi_clinical.AgeAtVisit(170) = NaN;
    genfi_clinical.EYO(170) = NaN;
    genfi_clinical.DiseaseDuration{170} = 'NA';

    clin_stages_table = readtable('/data/SD_FTD/jill/MCM_files/clin_stage_table_sensitivity.txt', 'ReadVariableNames',false);
    
end

clin_stages_table.Properties.VariableNames = {'ID', 'clin_stage', 'clin_subtype'};

clin_subtype_names = ["Control", "C9orf72", "GRN", "MAPT"];

%extract cTI pseudotime, subtypes, feature contributions and create tables for each folder (if cTI run for each factor)

clear pseudotime_column_names subtype_column_names N_trajs
folder_counter = 1;
for folder = 1:length(cti_folders)
    if cti_folders(folder).isdir

        % extract cTI pseudotimes and subtypes
        cti_file = dir([path filesep cti_folders(folder).name filesep 'cTI_IDs*.txt']);
        cti_data = readtable([path filesep cti_folders(folder).name filesep cti_file.name]);
        pseudotime_name = string(['cTI_pseudotime_' cti_folders(folder).name]);
        
        extract cTI feature contributions
        features_contributions_file = dir([path filesep cti_folders(folder).name filesep 'cTI_features_contributions*.txt']);
        features_contributions_data = load([path filesep cti_folders(folder).name filesep features_contributions_file.name]);
        features_name = string([cti_folders(folder).name]);
        features_contributions_data = table(features_contributions_data, 'VariableNames', features_name);
              
        %get subtype names
        N_cluster_rows = size(cti_data,2) - 2;
        clear subtype_names
        for cluster = 1:N_cluster_rows
            subtype_names(cluster,:) = string(['cTI_subtypes_' cti_folders(folder).name '_' num2str(cluster)]);
        end
        %get number of subtypes
        N_trajs(:,folder) = N_cluster_rows;
        
        % get indices of subjects in each subtype
        input_file = dir([path filesep cti_folders(folder).name filesep 'Input_data*.mat']);
        load([path filesep cti_folders(folder).name filesep input_file.name], 'cTI_pseudotemporal_paths');
        subtype_inds = zeros(size(cti_data,1), length(cTI_pseudotemporal_paths));
        for i = 1:length(cTI_pseudotemporal_paths)
            subtype_inds(1:length(cTI_pseudotemporal_paths(i).points), i) = cTI_pseudotemporal_paths(i).points;   
        end
        subtype_inds_data = array2table(subtype_inds, 'VariableNames', subtype_names');
        % name columns (ids, cTI stages and subtypes)
        cti_data.Properties.VariableNames = ["ID", pseudotime_name, subtype_names'];
        pseudotime_column_names(folder_counter,:) = pseudotime_name;
        subtype_column_names(folder_counter,:) = subtype_names(1,1);
        features_column_names(folder_counter,:) = features_name;
        
        %check if cTI ran correctly (if pseudotimes have any Inf values, something went wrong!)
        inf_ind = find(cti_data.(pseudotime_name) == Inf);
        if ~(isempty(inf_ind))
            disp([char(pseudotime_name) ' has Infinity values in cTI pseudotime'])
        end
        
        %combine data (cTI stages and subtypes) into table, create feature contributions table
        if folder_counter == 1
            cti_table = cti_data;
            subtype_inds_table = subtype_inds_data;
            subtypes_all_possible = cti_data(:, 3:end);
            features_contributions_table = features_contributions_data;
            N_features = size(features_contributions_table, 1);
        else
            cti_table = [cti_table cti_data(:,2:end)];
            subtype_inds_table = [subtype_inds_table subtype_inds_data];
            subtypes_all_possible = [subtypes_all_possible cti_data(:, 3:end)];
            %Create new table with row of empty values
            features_contributions_data{size(features_contributions_data,1)+1:N_features,:} = NaN;             
            features_contributions_table = [features_contributions_table features_contributions_data];

        end
        folder_counter = folder_counter + 1;
    end
end

% if cti only run once (ie no separate folders in path)
if ~(exist('pseudotime_column_names', 'var'))
    cti_file = dir([path filesep 'cTI_IDs*.txt']);
    cti_table = readtable([path filesep cti_file.name]);
    
    %get cTI feature contributions
    features_contributions_file = dir([path filesep 'cTI_features_contributions*.txt']);
    features_contributions_data = load([path filesep features_contributions_file.name]);
    features_name = "feature_contributions";
    features_contributions_table = table(features_contributions_data, 'VariableNames', features_name);
    
    %get subtype names
    N_cluster_rows = size(cti_table,2) - 2;
    clear subtype_names
    for cluster = 1:N_cluster_rows
        subtype_names(cluster,:) = string(['cTI_subtypes_' num2str(cluster)]);
    end
    %get number of subtypes
    %N_trajs = length(unique(cti_data(:,3)));

    cti_table.Properties.VariableNames = ["ID", "cTI_pseudotime", subtype_names'];
    pseudotime_column_names = "cTI_pseudotime";
    subtype_column_names = subtype_names(1,1);
    
    %check if pseudotimes have any Inf values
    inf_ind = find(cti_table.cTI_pseudotime == Inf);
    if ~(isempty(inf_ind))
        disp 'cTI_pseudotime has Infinity values in cTI pseudotime'
    end
end

% add clinical/EYO stages/subtypes to table
cti_table = [cti_table clin_stages_table(:,2:end)];

%add relevant columns from genfi_clinical table
cti_table.EYO = genfi_clinical.EYO;
duration_string = genfi_clinical(:, "DiseaseDuration");
duration_double = str2double(table2array(duration_string));
duration_double_table = array2table(duration_double, 'VariableNames', "DiseaseDuration");
clinical_strings = genfi_clinical(:, [47 58 70 72 74 76 78:84]);
clinical_variable_names = clinical_strings.Properties.VariableNames;
clinical_double = str2double(table2array(clinical_strings));
clinical_double_table = array2table(clinical_double, 'VariableNames', clinical_variable_names);
cti_table = [cti_table duration_double_table clinical_double_table];

% save table
writetable(cti_table, [path filesep 'cti_table.txt'])
% writetable(subtype_inds_table, [path filesep 'subtype_indices_table.txt'])

% add gene_status to cti_table
cti_table.Gene_Status = gene_status;

% get table of just carriers
cti_table_carriers = cti_table(cti_table.Gene_Status ~= 1, :);

% get matrix of pseudotime results
pseudotime_table = cti_table(:, pseudotime_column_names);
pseudotime = table2array(pseudotime_table);
pseudotime_names = pseudotime_table.Properties.VariableNames;

%% %%%%%%%%%%%%%%% Correlations 

%%% correlations between cTI and clinical stages/subtypes
[corr_stages, p_stages] = corr(pseudotime, cti_table.clin_stage, 'Type','Spearman');

%convert results to table
corr_stages_table = array2table(corr_stages, 'RowNames', pseudotime_names, 'VariableNames', "clinical_stage");
p_stages_table = array2table(p_stages, 'RowNames', pseudotime_names, 'VariableNames', "clinical_stage");


%%% correlations between cTI and clinical/neuropsych tests
[corr_clinical, p_clinical] = corr(pseudotime, clinical_double, 'Rows','pairwise');

% convert results to table
corr_clinical_table = array2table(corr_clinical, 'RowNames', pseudotime_names, 'VariableNames', clinical_variable_names);
p_clinical_table = array2table(p_clinical, 'RowNames', pseudotime_names, 'VariableNames', clinical_variable_names);

% find clinical r ranges
clear min_max_clinical_corr
for i = 1:length(pseudotime_names)
    correlations = table2array(corr_clinical_table(pseudotime_names(i),:));
    min_max_clinical_corr(i,:) = [min(abs(correlations)) max(abs(correlations))];
end
min_max_clinical_corr_table = array2table(min_max_clinical_corr, 'RowNames', pseudotime_names, 'VariableNames', ["Min" "Max"]);

%get only presymptomatic carriers
genfi_asymp = genfi_clinical(strcmp(genfi_clinical.GeneticStatus2, 'pos - asymp'), :);
pseudotime_asymp = pseudotime(strcmp(genfi_clinical.GeneticStatus2, 'pos - asymp'), :);
clinical_double_asymp = clinical_double(strcmp(genfi_clinical.GeneticStatus2, 'pos - asymp'), :);

%get only symptomatic carriers
genfi_symp = genfi_clinical(strcmp(genfi_clinical.GeneticStatus2, 'pos - symp'), :);
pseudotime_symp = pseudotime(strcmp(genfi_clinical.GeneticStatus2, 'pos - symp'), :);
clinical_double_symp = clinical_double(strcmp(genfi_clinical.GeneticStatus2, 'pos - symp'), :);

%exclude non-carriers
genfi_carriers = genfi_clinical(~strcmp(genfi_clinical.GeneticStatus2, 'neg'), :);
pseudotime_carriers = pseudotime(~strcmp(genfi_clinical.GeneticStatus2, 'neg'),:);
clinical_double_carriers = clinical_double(~strcmp(genfi_clinical.GeneticStatus2, 'neg'), :);

% correlations between cTI and clinical/neuropsych tests (presymp only)
[corr_clinical_asymp, p_clinical_asymp] = corr(pseudotime_asymp, clinical_double_asymp, 'Rows','pairwise');

% convert results to table
corr_clinical_asymp_table = array2table(corr_clinical_asymp, 'RowNames', pseudotime_names, 'VariableNames', clinical_variable_names);
p_clinical_asymp_table = array2table(p_clinical_asymp, 'RowNames', pseudotime_names, 'VariableNames', clinical_variable_names);

% correlations between cTI and clinical/neuropsych tests (symp only)
[corr_clinical_symp, p_clinical_symp] = corr(pseudotime_symp, clinical_double_symp, 'Rows','pairwise');

% convert results to table
corr_clinical_symp_table = array2table(corr_clinical_symp, 'RowNames', pseudotime_names, 'VariableNames', clinical_variable_names);
p_clinical_symp_table = array2table(p_clinical_symp, 'RowNames', pseudotime_names, 'VariableNames', clinical_variable_names);

% correlations between cTI and clinical/neuropsych tests (carriers only)
[corr_clinical_carriers, p_clinical_carriers] = corr(pseudotime_carriers, clinical_double_carriers, 'Rows','pairwise');

% convert results to table
corr_clinical_carriers_table = array2table(corr_clinical_carriers, 'RowNames', pseudotime_names, 'VariableNames', clinical_variable_names);
p_clinical_carriers_table = array2table(p_clinical_carriers, 'RowNames', pseudotime_names, 'VariableNames', clinical_variable_names);

%%% correlations between cTI and EYO
[corr_eyo, p_eyo] = corr(pseudotime, genfi_clinical.EYO, 'rows', 'complete');

% correlations between cTI and EYO, excluding non-carriers
[corr_eyo_carriers, p_eyo_carriers] = corr(pseudotime_carriers, genfi_carriers.EYO, 'rows', 'complete');

% correlations between cTI and EYO, only asymptomatic
[corr_eyo_asymp, p_eyo_asymp] = corr(pseudotime_asymp, genfi_asymp.EYO, 'rows', 'complete');

% correlations between cTI and EYO, only symptomatic
[corr_eyo_symp, p_eyo_symp] = corr(pseudotime_symp, genfi_symp.EYO, 'rows', 'complete');

% correlation between cTI and disease duration
[corr_duration, p_duration] = corr(pseudotime, duration_double, 'rows', 'complete');

%convert results to table
corr_eyo_table = table(corr_eyo, corr_eyo_carriers, corr_eyo_asymp, corr_eyo_symp, corr_duration, 'RowNames', pseudotime_names, 'VariableNames', ["EYO All", "EYO Carriers", "EYO Asymp", "EYO Symp", "Disease Duration"]);
p_eyo_table = table(p_eyo, p_eyo_carriers, p_eyo_asymp, p_eyo_symp, p_duration, 'RowNames', pseudotime_names, 'VariableNames', ["EYO All", "EYO Carriers", "EYO Asymp", "EYO Symp", "Disease Duration"]);
 

% save correlations to matlab file
save([path filesep 'correlations.mat'], 'corr_stages_table', 'p_stages_table', 'corr_clinical_table', ...
    'p_clinical_table', 'corr_eyo_table', 'p_eyo_table', 'min_max_clinical_corr_table', ...
    'corr_clinical_asymp_table', 'p_clinical_asymp_table', 'corr_clinical_carriers_table', ...
    'p_clinical_carriers_table', 'corr_clinical_symp_table', 'p_clinical_symp_table') 

% create table of correlations

stats_table = corr_clinical_carriers_table;
stats_table.EYO_carriers = corr_eyo_table.("EYO Carriers");
stats_table_transpose = rows2vars(stats_table);
stats_table_transpose.Properties.VariableNames{'OriginalVariableNames'} = 'Test';

p_table = p_clinical_table;
p_table.EYO_carriers = p_eyo_table.("EYO Carriers");
p_table_transpose = rows2vars(p_table);
p_table_transpose.Properties.VariableNames{'OriginalVariableNames'} = 'Test';

save([path filesep 'stats_tables.mat'], 'stats_table', 'p_table')


%% %%%%%%%%% Which data features predict better pseudotimes and subtrajetories?  

for i = 1:length(features_column_names)
    features_contributions = table2array(features_contributions_table(:,i));
    % All factors
    if i == 1
                
        %extract features for each modality
        clear features_mod_total
        for factor = 1:length(factor_names)
            j = factor - 1;
            if factor == 1
                ind = 1:N_regions;
            else
                ind = N_regions*j+1:N_regions*j+N_regions; %79:156, 157:234, 235:312, 313:390
            end
            features_mod = features_contributions(ind);
            
            %matrix of regions x modalities
            features_matrix(:,factor) = features_mod;
        end
        
        %sum across regions
        features_mod_total = sum(features_matrix);
        
        %sum across factors
        features_regions_total = sum(features_matrix, 2);
        
        % Modalities contribution across all features
        features_mod_table = table(factor_names', features_mod_total');
        features_mod_table_sorted = sortrows(features_mod_table,2,'descend');
        
	%save sorted table
        writetable(features_mod_table_sorted, [path filesep 'features_mod_table_sorted'])
        
        % Regional contribution
        features_regions_table = table(region_names, features_regions_total, features_matrix);
        features_regions_table_sorted = sortrows(features_regions_table,2,'descend');
        
        %save sorted table
        writetable(features_regions_table_sorted, [path filesep 'features_regions_table_sorted'])       
                  

%% %%%%%% Create Figure 1 for Paper

tiledlayout(2,2);
nexttile

%MMSE
col = 1;
%scatter(cti_table.MMSE, cti_table.(pseudotime_column_names(col)), 10, 'k', 'filled')
scatter(cti_table_carriers.MMSE, cti_table_carriers.(pseudotime_column_names(col)), 10, 'k', 'filled')
ylabel("cTI disease score")
yticks([0, 0.25, 0.50, 0.75, 1])
%xticks([0 10 20 30])
xlabel("MMSE")
hold on 
h = lsline;
h.Color = [0 0.4470 0.7410];
h.LineWidth = 2;
title('A                                                         ')

%CBI
nexttile
%scatter(cti_table.CBI_Total, cti_table.(pseudotime_column_names(col)), 10, 'k', 'filled')
scatter(cti_table_carriers.CBI_Total, cti_table_carriers.(pseudotime_column_names(col)), 10, 'k', 'filled')
ylabel("cTI disease score")
yticks([0, 0.25, 0.50, 0.75, 1])
%xticks([0 50 100])
xlabel("CBI")
hold on 
h = lsline;
h.Color = [0 0.4470 0.7410];
h.LineWidth = 2;
title('B                                                         ')

%by gene status (non-carrier, asymp, and symp)
nexttile
figure_name = ['gene_status_' char(pseudotime_column_names(col)) '.jpg'];
notBoxPlot(cti_table.(pseudotime_column_names(col)), cti_table.Gene_Status);
xticklabels({'non-carrier', 'presymptomatic', 'symptomatic'})
xtickangle(45)
ylabel("cTI disease score")
yticks([0, 0.25, 0.50, 0.75, 1])
title('C                                                         ')

%EYO
nexttile
scatter(genfi_carriers.EYO, pseudotime_carriers(:,1), 10, 'k', 'filled')
ylabel("cTI disease score")
yticks([0.00, 0.25, 0.50, 0.75, 1.00])
%xticks([-50 -25 0 25])
xlabel("EYO")
hold on 
h = lsline;
h.Color = [0 0.4470 0.7410];
h.LineWidth = 2;
title('D                                                         ')

