% Written by: Michael Lassi
% michael.lassi@santannapisa.it

if isfile(fullfile(PROCESSING_LIST.FILEPATH_OUTPUT{i},PROCESSING_LIST.FILENAME_MANUAL_CLEANING{i}))
   
    EEG = pop_loadset(PROCESSING_LIST.FILENAME_MANUAL_CLEANING{i}, PROCESSING_LIST.FILEPATH_OUTPUT{i});
    
    % Save original EEG to compute percent of signal removed
    old = EEG;
    
    command = '[EEG LASTCOM] = eeg_eegrej(EEG, eegplot2event(TMPREJ, -1));';
    eegplot(EEG.data, 'srate', EEG.srate,  'events', EEG.event, 'command', command,'winlength', 20);
    waitfor( findobj('parent', gcf, 'string', 'REJECT'), 'userdata');

    EEG.additional_info.removed_percent = (size(old.data,2) - size(EEG.data,2)) ./ size(old.data,2);

    % Here I save a first version of the data (manually rejected)
    pop_saveset_with_mkdir( EEG, FILENAME_MANUAL_CLEANING_PASS, FILEPATH_OUTPUT);

    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
end