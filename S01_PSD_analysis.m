%% PSD analysis script
% This script performs the PSD analysis on each subject and then plots the
% main results and statistics
% EEGLAB is a prerequisite to run the script

% EEG files should be already preprocessed

% Written by: Michael Lassi
% michael.lassi@santannapisa.it
%%
clear all
close all
clc
%% Start EEGlab
eeglab nogui

colors = [hex2rgb('#279084');hex2rgb('#FF7C70'); hex2rgb('#923A56')];
alphacolor = 0.3;
%% Setup
% EEG files should be saved with an identifier on their ID, that should
% match the csv file containing diagnostic label
BASE_FOLDER = ''; % Insert folder containing all info
EEG_FOLDER = [BASE_FOLDER,'']; % insert folder with preprocessed .set EEG files
file_list = getAllFiles(EEG_FOLDER,0,'.set'); % list all the eeg files
file_list = sort_nat(file_list);

% In our case the files are called 'EEG_PREVIEW_%ID%', so we extract the ID
% from there
subject_id = cellfun(@(x) str2num(x), extractBefore(extractAfter(file_list,'EEG_PREVIEW_'),FILE_ENDING));
OUTPUT = [BASE_FOLDER,'\PREVIEW_data\Processed\'];

% Sample path to one of teh eeg files to get info about recordings
% (sampling rate etc)
SAMPLE_EEG = ''; % insert sample EEG file
EEG = pop_loadset(SAMPLE_EEG);


%% PSD Computation

%%%%%%%%%%%%%%% SETTINGS FOR SPECTRA COMPUTATION %%%%%%%%%%%%%%%%%%%%%%%%%%
opts = [];
opts.spatial_filtering = '';
opts.bands = [1,4; 4 8; 8 13; 13 30];
opts.robust_metric = 0;
opts.full_range = [1 48];
opts.channel_norm = 0;
opts.PSD_WINDOW_LENGTH = 2; % seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for file = 1 : length(file_list)
    
    EEG = pop_loadset(file_list{file});
    EEG = pop_select( EEG, 'nochannel',{'ECG','EOG','64','MK','A1','A2'}); % unused channels
    disp('============')
    disp('Loaded file:')
    disp(file_list{file});
    disp('============')
    
    % Collect basic information
    PSD_WINDOW_LENGTH = 2;
    fs = EEG.srate;
    data = EEG.data';
    window_length = PSD_WINDOW_LENGTH * fs; % 2 s
    frequency_step = fs ./window_length;
    [psd(file,:,:), f] = pwelch(EEG.data', window_length,0,0:frequency_step:fs/2  , fs); % The result will be in muV^2 / Hz > convert to dB?
    
    % Extract relative power metrics from the psd of each channel
    opts.mode = 'normpower';
    [normpower_metrics(file,:),labels, channel_normpower(:,:,file)] = get_psd_metrics(squeeze(psd(file,:,:)),f,EEG,opts);
end
chanlocs = EEG.chanlocs;

%% Normality test
% Load the file with diagnostic information
GROUPS_FILE = [BASE_FOLDER, '']; % insert path to file here
group = get_group_label(GROUPS_FILE, subject_id);
%% Average power by group
% First, we check for normality by means of Shapiro-Wilk test.
% Then, we use ANOVA if the sample of the majority of variables are normal,
% otherwise we use kruskalwallis

% We only retain metrics we are interested in (average power by bands, by
% ROI)
labels_of_interest = [1:4, 13:36];
normpower_metrics = normpower_metrics(:,labels_of_interest);
labels = labels(labels_of_interest);

for i =1 : length(labels)
    area_name = labels{i};
    loi = find(strcmp(labels,area_name));
    
    metric_db = 10*log10(normpower_metrics(:,loi));
    % normality test
    [~,p_normality(i)] = swtest(metric_db);
end

if sum(p_normality > 0.05) >= length(p_normality)/2
    test_type = 'parametric';
else
    test_type = 'nonparametric';
end

disp(['Using ', test_type,' tests, given the results of the Shapiro-Wilk test']);
%% PSD metrics by group
% A label is HC
% B label is SCD
% C label is MCI

n_groups = unique(group);
SIGNIFICANCE_THRESHOLD = 0.05;
for i =1 : length(labels)
    
    area_name = labels{i};
    loi = find(strcmp(labels,area_name));
    
    % we set the dataset for statistics
    data_mat = prepare_unpaired_dataset(normpower_metrics(:,loi),group);
    % we convert the metrics to dB
    metric_db = 10*log10(normpower_metrics(:,loi));
    
    if strcmp(test_type, 'parametric')
        [p_groups(i),table,stats_groups] = anova1(metric_db,group,'off');
        f_groups(i) = table{2,5};
        [comparisons] = multcompare(stats_groups, SIGNIFICANCE_THRESHOLD, 'off','bonferroni');
        
    elseif strcmp(test_type,'nonparametric')
        [p_groups(i),~,stats_groups] = kruskalwallis(metric_db, group, 'off');
        [comparisons] = multcompare(stats_groups, SIGNIFICANCE_THRESHOLD, 'off','bonferroni');
    end
    eta_square(i) = get_eta_square(table);
    
    % Already Bonferroni corrected by multcompares
    idx_h = find(strcmp(stats_groups.gnames,'A'));
    idx_scd = find(strcmp(stats_groups.gnames,'B'));
    idx_mci = find(strcmp(stats_groups.gnames,'C'));
    
    [p_hscd(i)] = comparisons((comparisons(:,1) == idx_h & comparisons(:,2) == idx_scd) | (comparisons(:,2) == idx_h & comparisons(:,1) == idx_scd),6);
    [p_hmci(i)] = comparisons((comparisons(:,1) == idx_h & comparisons(:,2) == idx_mci) | (comparisons(:,2) == idx_h & comparisons(:,1) == idx_mci),6);
    [p_scdmci(i)] = comparisons((comparisons(:,1) == idx_scd & comparisons(:,2) == idx_mci) | (comparisons(:,2) == idx_scd & comparisons(:,1) == idx_mci),6);
    
    
    data_mat = 10*log10(data_mat);

    figure
    violinplot(data_mat);
    means(i,:) = mean(data_mat, 'omitnan');
    stds(i,:) = std(data_mat,[],'omitnan');
    d_hscd = get_cohens_d(data_mat(:,1), data_mat(:,2));
    d_hmci = get_cohens_d(data_mat(:,1), data_mat(:,3));
    d_scdmci = get_cohens_d(data_mat(:,2), data_mat(:,2));
    
    ylabel('Relative Power (dB)')
    xticklabels({'H','SCD','MCI'});
    xlabel('Group')
    
    if p_hmci(i) < SIGNIFICANCE_THRESHOLD
        sigstar([1 3], p_hmci(i))
    end
    if p_hscd(i) < SIGNIFICANCE_THRESHOLD
        sigstar([1 2], p_hscd(i))
    end
    
    if p_scdmci(i)  < SIGNIFICANCE_THRESHOLD
        sigstar([2 3], p_scdmci(i))
    end
    title([area_name],'Interpreter','none')
end
%% PSD plot by group
% Make plot of the psd spectrum by group with confidence interval
fig = figure;

options.x_axis =f;
options.error = 'c95';
options.handle = fig;

ALL_POWER = f >=1 & f <= 48;

options.line_width = 1.5;
options.alpha = 0.4;
k = 1;
grouping = unique(group);
for i =1 :length(grouping)
    
    options.color_area = colors(k,:);
    options.color_line = colors(k,:);
    power_norm = psd./sum(psd(:,ALL_POWER,:),2);
    
    plot_areaerrorbar(10*log10(mean(power_norm(strcmp(group, grouping{i}),:,:),3)), options);
    hold on
    xlim([1 30])
    k= k+1;
    xline(4,'k--');
    xline(8,'k--');
    xline(13,'k--');
    
end

set(gca,'box','off')
% adjust the legend
h3 = plot(NaN,NaN, 'Color',colors(1,:));
h1 = plot(NaN,NaN, 'Color', colors(2,:));
h2 = plot(NaN,NaN, 'Color', colors(3,:));
hl = legend([h3,h1,h2],{'H','SCD','MCI'});
legend boxoff
xlabel('Frequency (Hz)')
ylabel('Relative Power (dB)')
title('Average Relative Power')
%% Topographies by group
power_idx = 4; % the power (delta, theta, alpha, beta) to plot
power_name = {'delta','theta','alpha','beta'};

grouping = {'H','SCD','MCI'};
cmap = 'parula';
set_graphics_paper(3,2)

figure
groups = unique(group);
ngroups = length(groups);

for i = 1: ngroups
    subjects = strcmp(group,groups{i});
    power_to_plot(:,i)= squeeze(mean(channel_normpower(power_idx,:,subjects),3));
    absmin = min(power_to_plot(:));
    absmax = max(power_to_plot(:));
end

for i = 1 : ngroups
    subplot(1,ngroups,i)
    
    topoplot(power_to_plot(:,i), EEG.chanlocs,'maplimits',[absmin,absmax]);
    colormap(cmap)
    title(grouping{i})
    side_colorbar(1,ngroups, absmin, absmax, cmap);
    
end

