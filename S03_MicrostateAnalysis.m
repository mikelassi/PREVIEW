%%  Microstate analysis
% This script relies on EEGLAB and on Thomas Koenig Microstate toolbox to
% perform microstate extraction and averaging. For the backfitting, it uses
% the Microstate EEGlab toolbox

% Written by: Michael Lassi
% michael.lassi@santannapisa.it
clear all
close all
clc

eeglab nogui
%%
BASE_FOLDER = ''; % Insert folder with the entire project

OUTPUT = [BASE_FOLDER,'']; % insert folder where to save microstates data
FOLDER_PREP = [BASE_FOLDER,'']; % Insert folder with preprocessed eegs
files = getAllFiles(FOLDER_PREP,0,'.set');
files = sort_nat(files);


[~,filenames] = fileparts(files);
splitted = split(filenames,'_');
subject_id = str2double(splitted(:,4));
%% Select only indices of the single groups
GROUP_FILE = [BASE_FOLDER, '']; %insert path to file with diagnostic groups
[true_group] = get_group_label(GROUP_FILE,subject_id);
cond2diag = containers.Map({'A','B','C','NA'},{'H','SCD','MCI','NA'});

%% Load all the subjects in EEGLAB
sub_eeg = [];

for i = 1 : length(files)
    EEG = pop_loadset(files{i});
    
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG,EEG,0);
    ALLEEG = eeg_store(ALLEEG, EEG, CURRENTSET);
    sub_eeg = [sub_eeg CURRENTSET];
end

%% Now perform the clustering
PLOTTING = 1;
ClustPars = struct('MinClasses',4,'MaxClasses',4, 'GFPPeaks',true,...
    'IgnorePolarity',true,'MaxMaps',3e4,'Restarts',50,'UseAAHC',false);
chanlocs = EEG.chanlocs;

for i = 1 : length(files)
    tmpEEG = eeg_retrieve(ALLEEG, sub_eeg(i));
    tmpEEG = pop_FindMSTemplates(tmpEEG,ClustPars);
    
    if PLOTTING
        figure;plot_topo_set(tmpEEG.msinfo.MSMaps(4).Maps', chanlocs);
        sgtitle(escape_special_characters(['Subject ', num2str(sub_eeg(i))]))
        colormap jet
        save_image_hd(['subject_', num2str(sub_eeg(i)),'_microstatesKoenig'],'on', IMAGE_OUT);
    end
    close all
    
    ALLEEG = eeg_store(ALLEEG,tmpEEG,sub_eeg(i));
end
%% Write the group inside the EEG struct
for i = 1 : length(files)
    tmpid = strsplit(ALLEEG(i).setname,' ');
    subject_id(i) = str2num(tmpid{1});
    ALLEEG(i).subject = subject_id(i);
    
end

for i =1 : length(files)
    ALLEEG(i).group = cond2diag(true_group{i});
end

save('microstateKoenig.mat','-v7.3')

% load('microstateKoenig.mat')

%% Compute the grand mean and sort it manually
EEG = pop_CombMSTemplates(ALLEEG,sub_eeg, 0,0,'GrandMean');

[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'gui','off');
[ALLEEG, EEG] = pop_ShowIndMSMaps(EEG,4,1,ALLEEG);

GrandMeanIdx = CURRENTSET;

[ALLEEG, EEG] = pop_ShowIndMSMaps(EEG,4,1,ALLEEG);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);

%% SORT ALL MAPS BASED ON GRAND MEAN
ALLEEG = pop_SortMSTemplates(ALLEEG, sub_eeg, 0 , GrandMeanIdx);
ALLEEG = eeg_store(ALLEEG, EEG, CURRENTSET);
%%
healthy_idx = find(strcmp({ALLEEG.group},'H'));
EEG = pop_CombMSTemplates(ALLEEG,healthy_idx, 0,0,'MeanHealthy');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'gui','off');
healthy_mean_idx = CURRENTSET;
[ALLEEG, EEG] = pop_ShowIndMSMaps(EEG,4,1,ALLEEG);

%%
scd_idx = find(strcmp({ALLEEG.group},'SCD'));
EEG = pop_CombMSTemplates(ALLEEG,scd_idx, 0,0,'MeanSCD');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 5, 'gui','off');
scd_mean_idx = CURRENTSET;
[ALLEEG, EEG] = pop_ShowIndMSMaps(EEG,4,1,ALLEEG);

%% 
mci_idx = find(strcmp({ALLEEG.group},'MCI'));
EEG = pop_CombMSTemplates(ALLEEG,mci_idx, 0,0,'MeanMCI');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 5, 'gui','off');
mci_mean_idx = CURRENTSET;
[ALLEEG, EEG] = pop_ShowIndMSMaps(EEG,4,1,ALLEEG);

%% Extract maps and sorted maps and compare with the other method
for i = 1 : length(sub_eeg)
   eeg_analysis = ALLEEG(sub_eeg(i));
   F_u{i} = eeg_analysis.msinfo.MSMaps(4).Maps';
end
save('sorted_microstates.mat','F_u','subject_id');
%% Save grand means

grand_mean = ALLEEG(GrandMeanIdx).msinfo.MSMaps(4).Maps';
healthy_mean =  ALLEEG(healthy_mean_idx).msinfo.MSMaps(4).Maps';
scd_mean = ALLEEG(scd_mean_idx).msinfo.MSMaps(4).Maps';
mci_mean = ALLEEG(mci_mean_idx).msinfo.MSMaps(4).Maps';
save('sorted_means.mat', 'grand_mean', 'healthy_mean','scd_mean','mci_mean');

%% Topographical analysis
SAMPLE_EEG = ''; % Insert one preprocessed eeg to collect general info about recordings
EEG = pop_loadset(SAMPLE_EEG);
load('sorted_microstates_koenig.mat'); 
OUTPUT_FILE = fullfile(INPUT,'\topo_metrics_koenig.mat');
load('sorted_means.mat')
colors = [0.967797559291991,0.441274560091574,0.535810315505870; 0.504901784953007,0.590911923121528,0.958465725212856];
alphacolor = 0.3;
OUTPUT = '';% path to store the results

GROUP_FILE = ''; % File with diagnostic label
groups = get_group_label(GROUP_FILE,subject_id);

%% Grouping microstates by diagnosis
F_group = {{},{},{}};
for i =1 : length(F_u)
    if strcmp(groups{i}, 'A')
        F_group{1} = [F_group{1} F_u{i}];
    elseif strcmp(groups{i},'B')
        F_group{2} = [F_group{2} F_u{i}];
    elseif strcmp(groups{i},'C')
        F_group{3} = [F_group{3} F_u{i}];
    end
end

%% Sorting microstates
% Now I have the subject specific assignment for each group. I need to sort
% by that and then perform TANOVA
for i = 1 : length(F_group)
    for k =1 : size(F_group{i}{1},2)
        grouped_microstates{i}{k} = cell2mat(cellfun(@(x) x(:,k), F_group{i},'UniformOutput',false));
    end
end
%% TANOVA

grouping = [zeros(size(grouped_microstates{1}{1},2),1);ones(size(grouped_microstates{2}{1},2),1);2*ones(size(grouped_microstates{3}{1},2),1)];

for k =1 : size(grouped_microstates{1},2)
    [p_tan_groups(k), diss_tan_groups(k)] = tanova_multigroup([grouped_microstates{1}{k},grouped_microstates{2}{k},grouped_microstates{3}{k}]', grouping,10000,1);
end
for k =1 : size(grouped_microstates{1},2)
    [p_tan_scdmci(k), diss_tanova_scd_mci(k),diss_sh{k}] = tanova(grouped_microstates{2}{k}',grouped_microstates{3}{k}', 10000,1);
    
end

for k =1 : size(grouped_microstates{1},2)
    [p_tan_hscd(k), diss_tanova_h_scd(k),diss_sh{k}] = tanova(grouped_microstates{1}{k}',grouped_microstates{2}{k}', 10000,1);
    
end
for k =1 : size(grouped_microstates{1},2)
    [p_tan_hmci(k), diss_tanova_h_mci(k),diss_sh{k}] = tanova(grouped_microstates{1}{k}',grouped_microstates{3}{k}', 10000,1);
end

p_scdmci_bh = fwer_bonf(p_tan_scdmci,0.05)';
p_hmci_bh = fwer_bonf(p_tan_hmci,0.05)';
p_hscd_bh = fwer_bonf(p_tan_hscd,0.05)';


%% Extraction of time features from microstate topographies
% Extracts all the microstate time features from the EEG signals of the
% subjects. It assumes microstate grand mean was already obtained and it is
% used to backfit into single subjects signals.

eeglab nogui
%% Define path to relevant files
SMOOTHING_TIME = 0.03; %s
DATA_PATH = ''; % Path to were single subjects microstate are stored
IMAGE_PATH = '';
FOLDER_PREP =[DATA_PATH, '']; % path to eeg files
files = getAllFiles(FOLDER_PREP,0,'.set');

FILE_ENDING = '_final'; % final preprocessed files
files_final = files(contains(files,FILE_ENDING));
files_final = sort_nat(files_final);
subject_id = str2double(extractBefore(extractAfter(files_final,'EEG_PREVIEW_'),FILE_ENDING));

OUTPUT = [DATA_PATH,'']; % path to store results
GROUPS_FILE = [DATA_PATH, '']; % path with diagnostic groups
% This contains the sorted microstates for each subject
FILE_NAMING = '';
MIC_FILE = ['sorted_microstates_',FILE_NAMING,'.mat'];
% This contains the sorted microstates means
MEAN_MIC_FILE = 'sorted_means.mat';

%% Loading
% This is only needed for chanlocs
EEG_SAMPLE =[DATA_PATH,'']; % sample eeg file for channel locations
EEG = pop_loadset(EEG_SAMPLE);

true_group = get_group_label(GROUPS_FILE,subject_id);

load(MIC_FILE);
load(MEAN_MIC_FILE);
%% We need to do the smoothing with the MicroSmooth function, so we need to reload all the eegs
data = [];
opts = [];
opts.minTime = round(SMOOTHING_TIME *EEG.srate);
opts.polarity = 0;
%%
data = cell(length(files_final),1);
gfp = cell(length(files_final),1);
sequence = cell(length(files_final),1);

for i = 1 : length(files_final)
    EEG = pop_loadset(files_final{i});
    % remove A1 and A2 channels
    chan_idx = find(strcmp({EEG.chanlocs.labels},'A1'));
    chan_idx = [chan_idx find(strcmp({EEG.chanlocs.labels},'A2'))];
    
    EEG = pop_select(EEG, 'nochannel',chan_idx);
    data{i} =  EEG.data;
    gfp{i} = get_gfp(data{i});
    sequence{i}  = MicroSmooth(data{i}, grand_mean,'reject segments',opts);
    gev(i) = get_gev(sequence{i}, grand_mean, data{i});
end
% First checkpoint: save smoothed sequences
save(fullfile(OUTPUT,['MicroSmoothedSequence_', FILE_NAMING,'_',num2str(opts.minTime),'.mat']),'sequence', 'gfp','subject_id');

%% Reload smoothed sequences if you want to start back form here

% load(fullfile(OUTPUT,['MicroSmoothedSequence_15.mat']))
%% Extract  microstates statistics (appearance, transition probabilities...)

for i = 1 : length(sequence)
    disp(['Processing : ', num2str(i)])
    seq_duration = length(sequence{i}) ./ EEG.srate;
    sequence_reduced = reduce_sequence(sequence{i});
    appearances(i,:) = get_appearance(sequence_reduced);
    average_duration(i,:) = get_duration(sequence{i},EEG.srate);
    transition_probabilities(:,:,i) = get_transition_probabilities(sequence{i});
    transition_probabilities_reduced(:,:,i) = get_transition_probabilities(sequence_reduced);
    tmp = transition_probabilities_reduced(:,:,i);
    tpr(i,:) = tmp(:);
    tmp = transition_probabilities(:,:,i);
    tp(i,:) = tmp(:);
    lzzip(i,1) =get_zlib_complexity(sequence_reduced,2);
    frequency(i,:) = get_occurrence(sequence_reduced,seq_duration);
    hurst(i,1) = get_hurst_exponent(sequence{i});
end

save([DATA_PATH,'PREVIEW_data\Processed\micro_statistics_single.mat'],'lzzip','hurst', 'frequency','tpr', 'tp','appearances','average_duration', 'sequence');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% load([DATA_PATH,'PREVIEW_data\Processed\micro_statistics_singlemat']);
%%

mean_micro_duration = mean(average_duration,2);

%% Plots
colors = [hex2rgb('#279084');hex2rgb('#FF7C70'); hex2rgb('#923A56')];
alphacolor = 0.3;

xs = [lzzip,  hurst, appearances*100,frequency, average_duration*1000, mean_micro_duration*1000];
labels = {'LZ complexity','Hurst Exponent',...
     'Coverage (A) (%)','Coverage (B) (%)','Coverage (C) (%)','Coverage (D) (%)',...
   'Frequency (A) ','Frequency (B) ','Frequency (C) ', 'Frequency (D)', 'Duration (A) (ms)'...
    'Duration (B) (ms)','Duration (C) (ms)','Duration (D) (ms)','Mean Duration (ms)'};

%%
SIGNIFICANCE_THRESHOLD = 0.05;
% if the majority of the distributions is normal (sw test, apply
% parametric tests, otherwise use non-parametric ones);
for i =1 : length(labels)
    area_name = labels{i};
    loi = find(strcmp(labels,area_name));
    
    metric_db = xs(:,loi);
    % normality test
    [~,p_normality(i)] = swtest(metric_db);
end

if sum(p_normality > 0.05) > length(p_normality)/2
    test_type = 'parametric';
else
    test_type = 'nonparametric';
end

disp(['Using ', test_type,' tests, given the results of the Shapiro-Wilk test']);

%%
for i = 1 :size(xs,2)
    x = xs(:,i);%appearances * 100;%average_duration * 1000;
    set_graphics_paper(2,2)
    
    area_name = labels{i};
    loi = find(strcmp(labels,area_name));
    metric_db = xs(:,loi);
    
    if strcmp(test_type, 'parametric')
        [p_groups(i),~,stats_groups] = anova1(metric_db,true_group,'off');
        [comparisons] = multcompare(stats_groups, SIGNIFICANCE_THRESHOLD, 'off','bonferroni');
        
    elseif strcmp(test_type,'nonparametric')
        [p_groups(i),table,stats_groups] = kruskalwallis(metric_db, true_group, 'off');
        h_stat(i) = table{2,5};
        [comparisons] = multcompare(stats_groups, SIGNIFICANCE_THRESHOLD, 'off','bonferroni');
    end
    % Already corrected by multcompares
    idx_h = find(strcmp(stats_groups.gnames,'A'));
    idx_scd = find(strcmp(stats_groups.gnames,'B'));
    idx_mci = find(strcmp(stats_groups.gnames,'C'));
    
    [p_hscd(i)] = comparisons((comparisons(:,1) == idx_h & comparisons(:,2) == idx_scd) | (comparisons(:,2) == idx_h & comparisons(:,1) == idx_scd),6);
    [p_hmci(i)] = comparisons((comparisons(:,1) == idx_h & comparisons(:,2) == idx_mci) | (comparisons(:,2) == idx_h & comparisons(:,1) == idx_mci),6);
    [p_scdmci(i)] = comparisons((comparisons(:,1) == idx_scd & comparisons(:,2) == idx_mci) | (comparisons(:,2) == idx_scd & comparisons(:,1) == idx_mci),6);
    
    d_hscd(i) = get_cohens_d(data_matrix(:,1), data_matrix(:,2));
    d_hmci(i) = get_cohens_d(data_matrix(:,1), data_matrix(:,3));
    d_scdmci(i) = get_cohens_d(data_matrix(:,2), data_matrix(:,3));

    
    figure
 
    data_matrix = prepare_unpaired_dataset(x, true_group);
    medians(i,:) = median(data_matrix,'omitnan');
    iqrs(i,:) = iqr(data_matrix);
    
    v = violinplot(data_matrix(:,1:3));
    v(1).ViolinColor=  colors(1,:);
    v(2).ViolinColor = colors(2,:);
    v(3).ViolinColor = colors(3,:);

    if p_scdmci(i) < SIGNIFICANCE_THRESHOLD
        sigstar([2 3],p_scdmci(i));
    end
    if p_hmci(i) < SIGNIFICANCE_THRESHOLD
        sigstar([1 3],p_hmci(i));
    end
    if p_hscd(i) < SIGNIFICANCE_THRESHOLD
        sigstar([1 2],p_hscd(i));
    end
    xticklabels({'H','SCD','MCI'})
    ylabel(labels{i})  
end
%% Store transition probabilities in a readable way
group_name=  {'H','SCD','MCI'};
groups = unique(true_group);
ngroups = length(groups);

tpr_mat = reshape(tpr,length(sequence),4,4);

tpr_mat_sorted = permute(tpr_mat,[2,3,1]);
for i = 1 : ngroups
    tpr_plot(:,:,i) = mean(tpr_mat_sorted(:,:, strcmp(true_group, groups{i})),3);
    iqr_plot(:,:,i) = iqr(tpr_mat_sorted(:,:,strcmp(true_group, groups{i})),3);
    tpr_grouped{i} = tpr_mat_sorted(:,:, strcmp(true_group, groups{i}));
    
end

for i = 1 : size(tpr_grouped{1},1)
    for j = 1 : size(tpr_grouped{1},2)
        p_tr(i,j) = ranksum(squeeze(tpr_grouped{2}(i,j,:)), squeeze(tpr_grouped{3}(i,j,:)));
    end   
end

p_tr_lin = reshape(p_tr, 1,[]);
p_corrected_tr = reshape(fwer_bonf(p_tr_lin,0.05),4,4);
%% Perform statistics on transition probabilities
SIGNIFICANCE_THRESHOLD = 0.05;
abs_min = min(tpr_plot(:));
abs_max = max(tpr_plot(:));

set_graphics_paper(3.5,1)
figure
% and plot them
h_sub = tight_subplot(1,ngroups ,[.1 .1],[.2 .2],[0.15 .15]);

for i = 1:ngroups
    axes(h_sub(i));
    
    imagesc(tpr_plot(:,:,i))
    colormap inferno
    caxis([abs_min abs_max]);
    xticks(1:4)
    yticks(1:4);
    
    xticklabels({'A','B','C','D'})
    yticklabels({'A','B','C','D'})
    ylabel(' t ' )
    xlabel(' t + 1 ')
    title(group_name{i})
end

colorbar_labeled('')
sgtitle('Transition Probabilities Matrices')

% Statistical testing with Bonferroni correction
for i =1 : size(tpr_mat_sorted,1)
    for  j =1 : size(tpr_mat_sorted,2)
        
        metric_db = tpr_mat(:,i,j);
        if strcmp(test_type, 'parametric')
            [p_groups_tpr(i,j),~,stats_groups] = anova1(metric_db,true_group,'off');
            [comparisons] = multcompare(stats_groups, SIGNIFICANCE_THRESHOLD, 'off');
            
        elseif strcmp(test_type,'nonparametric')
            [p_groups_tpr(i,j),table,stats_groups] = kruskalwallis(metric_db, true_group, 'off');
            h_groups_tpr(i,j) = table{2,5};
            [comparisons] = multcompare(stats_groups, SIGNIFICANCE_THRESHOLD, 'off','bonferroni');
        end
        % Already corrected by multcompares
        
        idx_h = find(strcmp(stats_groups.gnames,'A'));
        idx_scd = find(strcmp(stats_groups.gnames,'B'));
        idx_mci = find(strcmp(stats_groups.gnames,'C'));
        d_tpr_hscd(i,j) = get_cohens_d(metric_db(strcmp(groups{1}, true_group),:),...
            metric_db(strcmp(groups{2}, true_group),:));
        d_tpr_hmci(i,j) = get_cohens_d(metric_db(strcmp(groups{1}, true_group),:),...
            metric_db(strcmp(groups{3}, true_group),:));
         d_tpr_scdmci(i,j) = get_cohens_d(metric_db(strcmp(groups{2}, true_group),:),...
            metric_db(strcmp(groups{3}, true_group),:));
        
        [p_hscd_tpr(i,j)] = comparisons((comparisons(:,1) == idx_h & comparisons(:,2) == idx_scd) | (comparisons(:,2) == idx_h & comparisons(:,1) == idx_scd),6);
        [p_hmci_tpr(i,j)] = comparisons((comparisons(:,1) == idx_h & comparisons(:,2) == idx_mci) | (comparisons(:,2) == idx_h & comparisons(:,1) == idx_mci),6);
        [p_scdmci_tpr(i,j)] = comparisons((comparisons(:,1) == idx_scd & comparisons(:,2) == idx_mci) | (comparisons(:,2) == idx_scd & comparisons(:,1) == idx_mci),6);
        
             
        if p_scdmci_tpr(i,j) < 0.05 && p_scdmci_tpr(i,j) >= 0.01
            sig_scdmci{i,j} = '*';
        elseif p_scdmci_tpr(i,j) < 0.01 && p_scdmci_tpr(i,j) >= 0.001
            sig_scdmci{i,j} = '**';
        elseif p_scdmci_tpr(i,j) < 0.001
            sig_scdmci{i,j} = '***';
        else
            sig_scdmci{i,j} = ' ';
        end
        
        if p_hmci_tpr(i,j) < 0.05 && p_hmci_tpr(i,j) >= 0.01
            sig_hmci{i,j} = '*';
        elseif p_hmci_tpr(i,j) < 0.01 && p_hmci_tpr(i,j) >= 0.001
            sig_hmci{i,j} = '**';
        elseif p_hmci_tpr(i,j) < 0.001
            sig_hmci{i,j} = '***';
        else
            sig_hmci{i,j} = ' ';
        end
        
        if p_hscd_tpr(i,j) < 0.05 && p_hscd_tpr(i,j) >= 0.01
            sig_hscd{i,j} = '*';
        elseif p_hscd_tpr(i,j) < 0.01 && p_hscd_tpr(i,j) >= 0.001
            sig_hscd{i,j} = '**';
        elseif p_hscd_tpr(i,j) < 0.001
            sig_hscd{i,j} = '***';
        else
            sig_hscd{i,j} = ' ';
        end
    end
end
figure
N = size(tpr_mat,2);
x = repmat(1:N,N,1); % generate x-coordinates
y = x'; % generate y-coordinates

imagesc(tpr_plot(:,:,3)' - tpr_plot(:,:,2)')

colormap parula
colorbar
title('Transition Probability Differences (MCI - SCD)')
xticklabels({'A','B','C','D'})
yticklabels({'A','B','C','D'})
ylabel(' t ' )
xlabel(' t + 1 ')
text(x(:), y(:), sig_scdmci', 'HorizontalAlignment', 'Center','FontSize',20)

figure

imagesc(tpr_plot(:,:,2)' - tpr_plot(:,:,1)')
colormap parula
colorbar
title('Transition Probability Differences (SCD - H)')
xticklabels({'A','B','C','D'})
yticklabels({'A','B','C','D'})
ylabel(' t ' )
xlabel(' t + 1 ')
text(x(:), y(:), sig_hscd', 'HorizontalAlignment', 'Center','FontSize',20)

figure
imagesc(tpr_plot(:,:,3)' - tpr_plot(:,:,1)')
colormap parula
colorbar
title('Transition Probability Differences (MCI - H)')
xticklabels({'A','B','C','D'})
yticklabels({'A','B','C','D'})
ylabel(' t ' )
xlabel(' t + 1 ')
text(x(:), y(:), sig_hmci', 'HorizontalAlignment', 'Center','FontSize',20)

%% Plot sample gfps

idx = 79; % subjects 11
srate = EEG.srate;
time = [1 : length(sequence{idx})] * 1/srate;
figure
plot_micro_activation(sequence{idx},gfp{idx}, time)
xlim([2 4])
legend({'A','B','C','D'});
set(gca,'box','off')

