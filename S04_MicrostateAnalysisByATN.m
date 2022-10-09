%% Microstate analysis by ATN classification
% This script performs the microstate topographical extraction and analysis
% based on the groups identified by ATN classification. Specifically, we test A+ and A- subjects 
% and AD-pathology profile vs no AD pathology profile
% The needed toolobxes are the same as in S03

% Written by: Michael Lassi
% michael.lassi@santannapisa.it
clear all
close all
clc
eeglab nogui
%%
BASE_FOLDER = '';
IMAGE_OUT = '';
BASE_IMAGE = '';
OUTPUT = [BASE_FOLDER,'\PREVIEW_data\Processed\'];
FOLDER_PREP = [BASE_FOLDER,'\PREVIEW_data\Processed\EEG\EC\'];
files = getAllFiles(FOLDER_PREP,0,'.set');
files = sort_nat(files);

%
[~,filenames] = fileparts(files);
splitted = split(filenames,'_');
subject_id = str2double(splitted(:,4));

%% ABeta file
% 
abeta = readtable('');% insert file with abeta, p-tau and total-tau values

aplus = table2array(abeta(:,'ABeta')) < 670; 
tplus = table2array(abeta(:,'P_Tau')) > 60;
nplus = table2array(abeta(:,'Tau')) > 670;

diagnosis = table2array(abeta(:,'Diagnosi'));
% add the 
% Divide into abetamyloid-pathology vs not
abeta_pat = aplus;
% Divide into neurodegeneration en cours) vs not (A+(T/N +))
abeta_neurodeg = aplus & (tplus | nplus);


% Only A+
disp(['There are ',num2str(sum(strcmp(diagnosis,'SCD') & aplus)) ,' A+ SCD subjects ']);
disp(['There are ',num2str(sum(strcmp(diagnosis,'MCI') & aplus)) ,' A+ MCI subjects ']);
% AD-pathology
disp(['There are ',num2str(sum(strcmp(diagnosis,'SCD') & aplus & (tplus | nplus))) ,' AD-pathology SCD subjects ']);
disp(['There are ',num2str(sum(strcmp(diagnosis,'MCI') & aplus & (tplus | nplus))) ,' AD-pathology SCD subjects ']);

%% Load the single subject microstates computed in S03

load('microstate.mat')

%% Only keep subjects having the CSF measures

sub_eeg_abeta = table2array(abeta(:,'Numero_inOrdineDiRegistrazione_'));
abeta_flag = ismember(subject_id, sub_eeg_abeta);
ALLEEG(~abeta_flag) = [];

%% Compute the grand mean and sort it manually
sub_eeg =1:length(ALLEEG);
EEG = pop_CombMSTemplates(ALLEEG,sub_eeg, 0,0,'GrandMeanAbeta');

[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'gui','off');
[ALLEEG, EEG] = pop_ShowIndMSMaps(EEG,4,1,ALLEEG);

GrandMeanIdx = CURRENTSET;

[ALLEEG, EEG] = pop_ShowIndMSMaps(EEG,4,1,ALLEEG);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);

%% SORT ALL MAPS BASED ON GRAND MEAN
ALLEEG = pop_SortMSTemplates(ALLEEG, sub_eeg, 0 , GrandMeanIdx);
ALLEEG = eeg_store(ALLEEG, EEG, CURRENTSET);
%% Now divide by abeta vs non-abeta
abeta_idx = find(abeta_pat);
non_abeta_idx = find(~abeta_pat);

EEG = pop_CombMSTemplates(ALLEEG,abeta_idx, 0,0,'MeanAbeta');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'gui','off');
abeta_mean_idx = CURRENTSET;
[ALLEEG, EEG] = pop_ShowIndMSMaps(EEG,4,1,ALLEEG);

EEG = pop_CombMSTemplates(ALLEEG,non_abeta_idx, 0,0,'MeanNonAbeta');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'gui','off');
non_abeta_mean_idx = CURRENTSET;
[ALLEEG, EEG] = pop_ShowIndMSMaps(EEG,4,1,ALLEEG);

%% Degeneration vs non-degeneration
abeta_neurodeg_idx = find(abeta_neurodeg);
non_abeta_neurodeg_idx = find(~abeta_neurodeg);

EEG = pop_CombMSTemplates(ALLEEG,abeta_neurodeg_idx, 0,0,'MeanDegeneration');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'gui','off');
degeneration_mean_idx = CURRENTSET;
[ALLEEG, EEG] = pop_ShowIndMSMaps(EEG,4,1,ALLEEG);

EEG = pop_CombMSTemplates(ALLEEG,non_abeta_neurodeg_idx, 0,0,'MeanNonDegeneration');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'gui','off');
non_degeneration_mean_idx = CURRENTSET;
[ALLEEG, EEG] = pop_ShowIndMSMaps(EEG,4,1,ALLEEG);


%% Extract maps and sorted maps and compare with the other method
for i = 1 : length(sub_eeg)
   eeg_analysis = ALLEEG(sub_eeg(i));
   F_u{i} = eeg_analysis.msinfo.MSMaps(4).Maps';
end
save('sorted_microstates_abeta.mat','F_u','sub_eeg_abeta');
%% Save grand means

grand_mean = ALLEEG(GrandMeanIdx).msinfo.MSMaps(4).Maps';
abeta_mean =  ALLEEG(abeta_mean_idx).msinfo.MSMaps(4).Maps';
non_abeta_mean = ALLEEG(non_abeta_mean_idx).msinfo.MSMaps(4).Maps';
degeneration_mean = ALLEEG(degeneration_mean_idx).msinfo.MSMaps(4).Maps';
non_degeneration_mean = ALLEEG(non_degeneration_mean_idx).msinfo.MSMaps(4).Maps';

save('sorted_means_abeta.mat', 'grand_mean','abeta_mean',...
    'non_abeta_mean','degeneration_mean','non_degeneration_mean');


%%
BASE_FOLDER = ''; % general data folder
INPUT = [BASE_FOLDER,''];% folder containing the microstates files
SAMPLE_EEG  = ''; % sample folder to obtain eeg statistics
EEG = pop_loadset(SAMPLE_EEG);

colors = [0.967797559291991,0.441274560091574,0.535810315505870; 0.504901784953007,0.590911923121528,0.958465725212856];
alphacolor = 0.3;
OUTPUT = '';% output folder

%% load also the classic microstates
load('sorted_microstates.mat');
load('sorted_means.mat');

GROUP_FILE = 'D:\OneDrive - Scuola Superiore Sant''Anna\Data\PREVIEW_data\Raw\patients_groups.xlsx';
groups = get_group_label(GROUP_FILE, subject_id);

micro_healthy = F_u(strcmp(groups,'A'));

%% 
load('sorted_microstates_abeta.mat'); 
OUTPUT_FILE = fullfile(INPUT,'\topo_metrics_abeta.mat');
load('sorted_means_abeta.mat')

GROUP_FILE = 'abetapathology.xlsx';

groups = get_group_label(GROUP_FILE, sub_eeg_abeta);

group_names = unique(groups);
group_names = {'NonAbeta','Abeta'};
%% Divide microstates by groups
F_group = {{},{},{}};
for i =1 : length(F_u)
    for j =1 : length(group_names)
        if strcmp(groups{i}, group_names{j})
           F_group{j} = [F_group{j}, F_u{i}];
           continue
        end 
    end
end
%% Now perform the matching

for i = 1 : length(F_group)
    for k =1 : size(F_group{i}{1},2)
        grouped_microstates{i}{k} = cell2mat(cellfun(@(x) x(:,k), F_group{i},'UniformOutput',false));
    end
end
%% TANOVA

for k =1 : size(grouped_microstates{1},2)
    [p(k), diss(k),diss_sh{k}] = tanova(grouped_microstates{1}{k}',grouped_microstates{2}{k}', 10000,1); 
end
corrected_tanova = fwer_holmbonf(p,0.05)';