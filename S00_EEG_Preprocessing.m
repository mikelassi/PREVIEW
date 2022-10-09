%% Preprocessing script
% This performs PREP Pipeline + ICA removal of artifacts

% This script relies on EEGLAB with several add-on toolboxes to perform the
% preprocessing.
% specifically, the following toolboxes are needed: - iclabel, dipfit,
% binica

% Written by: Michael Lassi
% michael.lassi@santannapisa.it

%% Load subjects folder and EEG path
clear all
close all
clc

%% DEFINE BASE FOLDER AND PROCESSING FOLDER
BASE_FOLDER = ''; % Folder containing code, data all data and outputs location
CODE_FOLDER = [BASE_FOLDER,'']; % folder containing code

% Modify as needed
RAW_DATA_FOLDER = ''; % folder containing raw data files (edf files)
PROCESSED_DATA_FOLDER =''; % folder containg processed data as outputs

TMP_FOLDER = ''; % where to store tmp ica files (must not contain spaces)
%% Move to the working folder
cd(CODE_FOLDER);
addpath(genpath(CODE_FOLDER));


%% Initialize eeglab
% initialize eeglab structures
eeglab nogui;

%Force eeglab to load eegs as double instead of single variables (make rank
%computations more precise)
pop_editoptions( 'option_single', 0);

EEG_FOLDER = [BASE_FOLDER,RAW_DATA_FOLDER];
EEG_OUTPUT_FOLDER = [BASE_FOLDER, PROCESSED_DATA_FOLDER];

copy_folder_struct(EEG_FOLDER,EEG_OUTPUT_FOLDER);

SUBJECTS_FILES = getAllFiles(EEG_FOLDER,1,'.edf');

PARTIAL_FILESKIPPER = 1; % Select 1 if you want to start the analysis from where it was left

%% Initialize PROCESSING LIST STRUCTURE

[filepath, filename, extension] = fileparts(SUBJECTS_FILES);
new_filepath = strcat(EEG_OUTPUT_FOLDER,erase(strcat(filepath,'\'),EEG_FOLDER));

PROCESSING_LIST.SUBJECT_ID =cellstr( extractAfter(filename,'PREVIEW_'));
PROCESSING_LIST.FILE = cellstr(strcat(filepath,'\',filename,extension));
PROCESSING_LIST.FILEPATH_OUTPUT = cellstr(new_filepath);

PROCESSING_LIST.FILENAME_MANUAL_CLEANING = cellstr(strcat(filename,'_manual.set'));
PROCESSING_LIST.FILENAME_MANUAL_CLEANING_PASS = cellstr(strcat(filename,'_manual_pass.set'));

PROCESSING_LIST.FILENAME_ICA =cellstr(strcat(filename,'_ica.set'));
PROCESSING_LIST.FILENAME_PREPROCESSED =cellstr(strcat(filename,'_preprocessed.set'));
PROCESSING_LIST.FILENAME_FINAL = cellstr(strcat(filename,'_final.set'));

%% Data preparation, rereferencing and manual cleaning
for i = 1 : length(PROCESSING_LIST.FILE)
    FILEPATH_OUTPUT  = PROCESSING_LIST.FILEPATH_OUTPUT{i};
    FILENAME_MANUAL_CLEANING =  PROCESSING_LIST.FILENAME_MANUAL_CLEANING{i};
    
    if PARTIAL_FILESKIPPER && exist(fullfile(FILEPATH_OUTPUT,FILENAME_MANUAL_CLEANING),'file')
        continue;
    else
        disp('================================================================');
        disp(['Processing file: ',PROCESSING_LIST.FILENAME_MANUAL_CLEANING{i}]);
        disp('================================================================');
        EEGLAB_PREP_PIPELINE;
    end
end
% I would introduce a first manual cleaning here
%%  Manual cleaning
for i = 1 : length(PROCESSING_LIST.FILE) % Since you added this later, check for the presence of FILENAME FINAL AS WELL
    FILEPATH_OUTPUT  = PROCESSING_LIST.FILEPATH_OUTPUT{i};
    FILENAME_MANUAL_CLEANING_PASS =  PROCESSING_LIST.FILENAME_MANUAL_CLEANING_PASS{i};
    
    if PARTIAL_FILESKIPPER && (exist(fullfile(FILEPATH_OUTPUT,FILENAME_MANUAL_CLEANING_PASS),'file') ...
            || exist(fullfile(FILEPATH_OUTPUT,PROCESSING_LIST.FILENAME_FINAL{i})))
        continue;
    else
        EEGLAB_MANUAL_CLEANING_FIRST_PASS
    end
end
%% ICA
for i = 1 : length(PROCESSING_LIST.FILE)
    FILEPATH_OUTPUT  = PROCESSING_LIST.FILEPATH_OUTPUT{i};
    FILENAME_ICA =  PROCESSING_LIST.FILENAME_ICA{i};
    
    if PARTIAL_FILESKIPPER && exist(fullfile(FILEPATH_OUTPUT,FILENAME_ICA),'file')
        continue;
    else
        EEGLAB_CLEANING_AND_ICA
    end
end

%% ICA pruning
for i = 1 : length(PROCESSING_LIST.FILE)
    FILEPATH_OUTPUT  = PROCESSING_LIST.FILEPATH_OUTPUT{i};
    FILENAME_PREPROCESSED =  PROCESSING_LIST.FILENAME_PREPROCESSED{i};
    
    if PARTIAL_FILESKIPPER && exist(fullfile(FILEPATH_OUTPUT,FILENAME_PREPROCESSED),'file')
        continue;
    else
        EEGLAB_ICA_PRUNING_MANUAL
    end
end
%% FINAL MANUAL CLEANING OF RESIDUAL ARTIFACTS
for i = 1 : length(PROCESSING_LIST.FILE)
    FILEPATH_OUTPUT  = PROCESSING_LIST.FILEPATH_OUTPUT{i};
    FILENAME_FINAL =  PROCESSING_LIST.FILENAME_FINAL{i};
    
    if PARTIAL_FILESKIPPER && exist(fullfile(FILEPATH_OUTPUT,FILENAME_FINAL),'file')
        warning('File was already processed!');
        continue;
    else
        % We manually check for outlier eegs
        EEGLAB_MANUAL_CLEANING
    end
end
% END OF PREPROCESSING
