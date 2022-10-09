if isfile(fullfile(PROCESSING_LIST.FILEPATH_OUTPUT{i},PROCESSING_LIST.FILENAME_PREPROCESSED{i}))
    % Manually cleaning the residual dataset
    EEG = pop_loadset(PROCESSING_LIST.FILENAME_PREPROCESSED{i}, PROCESSING_LIST.FILEPATH_OUTPUT{i});
    command = '[EEG LASTCOM] = eeg_eegrej(EEG, eegplot2event(TMPREJ, -1));';
    eegplot(EEG.data, 'srate', EEG.srate,  'events', EEG.event, 'command', command,'winlength', 20);
    waitfor( findobj('parent', gcf, 'string', 'REJECT'), 'userdata');

    % Here I save a clean version of the data 
    pop_saveset_with_mkdir( EEG, FILENAME_FINAL, FILEPATH_OUTPUT);
    EEG = pop_delset( EEG, [1] );
    
end