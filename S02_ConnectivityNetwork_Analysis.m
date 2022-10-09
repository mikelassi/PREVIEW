%% Connectivity and network analysis
% This script performs the source reconstruction, connectivity
% estimation and extraction of  network metrics for each subject EEG data
% EEG should be provided as preprocessed .set files
% EEGLAB, fieldtrip and brain connectivity toolbox are needed to run the
% code

% Written by: Michael Lassi
% michael.lassi@santannapisa.it
%% Script to run the connectivity analysis on all patients
clear all
close all
clc
%% Source reconstruction and connectivity estimation
DATA_FOLDER = '';% Insert folder with preprocessed EEG files
OUTPUT_FOLDER =''; % insert folder to save results files
file_list = getAllFiles(DATA_FOLDER,1,'.set');
file_list = sort_nat(file_list);

copy_folder_struct(DATA_FOLDER,OUTPUT_FOLDER);

FILE_ENDING = 'final';
file_list(~contains(file_list,FILE_ENDING)) = [];
output = build_output_path(file_list, DATA_FOLDER,OUTPUT_FOLDER,'_icoh_aal.mat');
PARTIAL_FILESKIPPER = 1;

for i =1 : length(file_list)
    
    disp('=================================================================')
    disp(['Processing: ', file_list{i}]);
    disp(['File n. ', num2str(i),'/', num2str(length(file_list))]);
    disp('=================================================================')
    
    if exist(output{i},'file') && PARTIAL_FILESKIPPER
        warning([output{i},' already generated'])
        continue
    end
    % Perform the connectivity
    [sig_conn_bands, conn_pruned, cohspctrm_sh,parc,conn,bands,coh_bands,coh_bands_sh, sig_conn ]  = run_connectivity_coherence(file_list{i});
    save(output{i},'parc',...
        'coh_bands','sig_conn_bands','coh_bands_sh','sig_conn','bands','-v7.3');
    
    clear sig_conn sig_conn_bands conn_pruned cohspctrm_sh parc conn bands coh_bands coh_bands_sh
end

%% Convert the multiple connectivity files into a single one
CONN_FOLDER = ''; % insert connectivity folder (the output of previous step)
OUTPUT =''; % insert file to save the output to
file_list = getAllFiles(CONN_FOLDER,1,'.mat');

[~,filenames] = fileparts(file_list);
filter= contains(filenames,'aal');

file_list = file_list(filter);
filenames = filenames(filter);

file_list = sort_nat(file_list);
filenames = sort_nat(filenames);

splitted = split(filenames,'_');
subject_id = str2double(splitted(:,4));
textprogressbar('Loading Connectivity matrices: ')

for i = 1 : length(file_list)
    textprogressbar(i / length(file_list)*100)
    load(file_list{i});
    
    conn_matrix_complete(:,:,:,i) = coh_bands;
    
end

roi_info = parc.brainordinate;
save(OUTPUT,'conn_matrix_complete','bands','roi_info','subject_id','-v7.3')
textprogressbar('Completed!')
%% Estimate the network metrics (after pruning of the connectivity matrices)

PRUNING_MODE = 'connected';
CONN_FILE = 'D:\OneDrive - Scuola Superiore Sant''Anna\Data\PREVIEW_data\Processed\coherence.mat';
OUTPUT = ['D:\OneDrive - Scuola Superiore Sant''Anna\Data\PREVIEW_data\Processed\network_metrics_',PRUNING_MODE,'.mat'];
OUTPUT_ML = 'D:\OneDrive - Scuola Superiore Sant''Anna\Data\PREVIEW_data\\network_ml.mat';
load(CONN_FILE)
%% Select adjacency matrix
% 1st and second dimensions are adjacency matrices
% 3rd dimension is frequency band
% 4th dimension is subject

for idx = 1 : size(conn_matrix,4)
    disp(['Subject ', num2str(idx)])
    adj_mat = conn_matrix(:,:,:,idx);
    
    for band_idx = 1 : size(conn_matrix,3)
        disp(['Band ', num2str(band_idx)])
        
        adj_band = adj_mat(:,:,band_idx);
        adj_band = get_sparsified_matrix(adj_band, PRUNING_MODE, 0.3);
        [results(:, band_idx,idx), labels] = get_graph_measures(adj_band);
       
    end
end

save(OUTPUT, 'results', 'bands','labels','subject_id')

PRUNING_METHOD = 'connected';
GROUP_FILE = 'D:\OneDrive - Scuola Superiore Sant''Anna\Data\PREVIEW_data\Raw\patients_groups.xlsx';
CONNECTIVITY_FILE = ['D:\OneDrive - Scuola Superiore Sant''Anna\Data\PREVIEW_data\Processed\network_metrics_',PRUNING_METHOD,'.mat'];
IMAGE_OUT = 'D:\OneDrive - Scuola Superiore Sant''Anna\Images\PREVIEW\';
load(CONNECTIVITY_FILE);

groups = get_group_label(GROUP_FILE, subject_id);
% remove unused metrics
results(5,:,:) = [];
results(1,:,:) = [];

band_names = {'Delta'};
group_names = {'H','SCD','MCI'};
labels = {'Strength','Clustering Coefficient','Characteristic Path Length','Small-worldness'};

%% Testing normality of the variables
k = 1;
for i = 1 : length(labels)
    for j = 1 : length(band_names)
        metric_db = squeeze(results(i,j,:));
        % normality test
        [~,p_normality(k)] = swtest(metric_db);
        k = k + 1;
    end
end

if sum(p_normality > 0.05) >= length(p_normality)/2
    test_type = 'parametric';
else
    test_type = 'nonparametric';
end

disp(['Using ', test_type,' tests, given the results of the Shapiro-Wilk test']);
%% Plot and compute statistics
% Here I use a one-way ANOVA if data is normal, or a Kruskal-Wallis test if
% data is not normal
colors = [hex2rgb('#279084');hex2rgb('#FF7C70'); hex2rgb('#923A56')];
alphacolor = 0.3;

SIGNIFICANCE_THRESHOLD = 0.05;
figure
for i = 1 : length(labels)
    clear p
    
    for j = 1 : length(band_names)
        
        dataset = prepare_unpaired_dataset(squeeze(results(i,j,:)),groups);
        
        % Since connectivity displays normal overall behavior
        metric_db = squeeze(results(i,j,:));
        
        if strcmp(test_type, 'parametric')
            [p_groups(i,j),table,stats_groups] = anova1(metric_db,groups,'off');
            [comparisons] = multcompare(stats_groups, SIGNIFICANCE_THRESHOLD, 'off','bonferroni');
            
        elseif strcmp(test_type,'nonparametric')
            [p_groups(i,j),table,stats_groups] = kruskalwallis(metric_db, groups, 'off');
            [comparisons] = multcompare(stats_groups, SIGNIFICANCE_THRESHOLD, 'off','bonferroni');
        end
        f_groups(i,j) = table{2,5};
        eta_square(i,j) = get_eta_square(table);
        % Already Bonferroni corrected by multcompares
        idx_h = find(strcmp(stats_groups.gnames,'A'));
        idx_scd = find(strcmp(stats_groups.gnames,'B'));
        idx_mci = find(strcmp(stats_groups.gnames,'C'));
        
        [p_hscd(i,j)] = comparisons((comparisons(:,1) == idx_h & comparisons(:,2) == idx_scd) | (comparisons(:,2) == idx_h & comparisons(:,1) == idx_scd),6);
        [p_hmci(i,j)] = comparisons((comparisons(:,1) == idx_h & comparisons(:,2) == idx_mci) | (comparisons(:,2) == idx_h & comparisons(:,1) == idx_mci),6);
        [p_scdmci(i,j)] = comparisons((comparisons(:,1) == idx_scd & comparisons(:,2) == idx_mci) | (comparisons(:,2) == idx_scd & comparisons(:,1) == idx_mci),6);
        
        
        subplot(1,length(labels),i)
        v = violinplot(dataset);
        v(1).ViolinColor=  colors(1,:);
        v(2).ViolinColor = colors(2,:);
        v(3).ViolinColor = colors(3,:);
        medians(i,:) = median(dataset,'omitnan');
        iqrs(i,:) = iqr(dataset);
        d_hscd(i) = get_cohens_d(dataset(:,1), dataset(:,2));
        d_hmci(i) = get_cohens_d(dataset(:,1), dataset(:,3));
        d_scdmci(i) = get_cohens_d(dataset(:,2), dataset(:,3));
        
        xticklabels(group_names)
        if p_scdmci(i,j) <= 0.05
            sigstar([2 3], p_scdmci(i,j))
        end
        if p_hmci(i,j) <= 0.05
            sigstar([1 3], p_hmci(i,j))
        end
        if p_hscd(i,j) <= 0.05
            sigstar([1 2], p_hscd(i,j))
        end
        ylabel(escape_special_characters(labels{i}))   
    end   
end

