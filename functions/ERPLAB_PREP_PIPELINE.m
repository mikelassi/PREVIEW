% Written by: Michael Lassi
% michael.lassi@santannapisa.it

EEG.etc.eeglabvers = '2021.0'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = pop_fileio(PROCESSING_LIST.FILE{i}, 'dataformat','auto');

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',PROCESSING_LIST.SUBJECT_ID{i},'gui','off');

EEG = eeg_checkset( EEG );
EEG.data = double(EEG.data);

%% Removal of unused channels
EEG = pop_select( EEG, 'nochannel',{'ECG','EOG','64','MK','A1','A2','A','B'});
EEG = eeg_checkset( EEG );

%% Set channel locations to MNI coordinate file for BEM dipfit model
% Remove EEG from channel names
for j = 1 : length(EEG.chanlocs)
    channame = split(EEG.chanlocs(j).labels,' ');
    EEG.chanlocs(j).labels = channame{end};
end

EEG = pop_chanedit(EEG, 'lookup',[matlabroot,'\\toolbox\\eeglab\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc'],'convert','chancenter');
EEG = eeg_checkset( EEG );

%% Remove events (if present) and substitute with corresponding .txt file

% removal of 'empty' event
if ~isempty(EEG.event)
    EEG = pop_editeventvals(EEG,'delete',find(strcmp({EEG.event.type},'empty')));
end

% OLD Code for loading from event file (if no EDF+ available)
% [path,file] = fileparts(PROCESSING_LIST.FILE{i});
% EEG = load_events_from_file([file,'.txt'],path, EEG);



%% PREP pipeline
% this performs: 
% - drift removal 
% - line frequency removal (Cleanline)
% - robust average referencing
% - bad channel interpolation
EEG = pop_prepPipeline(EEG,struct('ignoreBoundaryEvents', true, 'lineFrequencies',...
    [50  100 150 200 250], 'reportMode', 'skipReport',  'consoleFID',1,...
    'cleanupReference', true, 'keepFiltered', true, 'removeInterpolatedChannels', false));


%% Here I save a first version of the data (manually rejected)
pop_saveset_with_mkdir( EEG, FILENAME_MANUAL_CLEANING, FILEPATH_OUTPUT);
EEG = pop_delset( EEG, [1] );

