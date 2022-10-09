% Written by: Michael Lassi
% michael.lassi@santannapisa.it

if isfile(fullfile(PROCESSING_LIST.FILEPATH_OUTPUT{i}, PROCESSING_LIST.FILENAME_MANUAL_CLEANING_PASS{i}))
    
    EEG = pop_loadset(PROCESSING_LIST.FILENAME_MANUAL_CLEANING_PASS{i}, PROCESSING_LIST.FILEPATH_OUTPUT{i}); 
    EEG.data = double(EEG.data);
    
    % To avoid binica complaining
    cd(TMP_FOLDER);
    %% ICA decomposition and component removal (manual)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% WARNING: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ica computation can extract at maximum m components, where m is the
    % rank of the EEG data. When removing and interpolating channels, we
    % are reducing eeg matrix rank to n_channels - k, where k is the number
    % of deleted channels. Binica should automatically detect this, but
    % sometimes it leads to weird and unstable results (for reference, see:
    % https://eeglab.org/tutorials/06_RejectArtifacts/RunICA.html#how-to-deal-with-corrupted-ica-decompositions)
    % hence, it is better to perform a pca before ica in those cases.
    % However, check that rank is computed correctly, as the function is
    % not super-solid in Matlab. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    eeg_rank = rank(EEG.data*EEG.data');
    if eeg_rank == length(EEG.chanlocs)
        EEG = pop_runica(EEG, 'icatype', 'binica', 'extended',1);
    else
        EEG = pop_runica(EEG, 'icatype', 'binica', 'pca', eeg_rank);
    end
    
    % Use ICLabel to classify the components
    EEG = eeg_checkset( EEG );
    EEG = pop_iclabel(EEG, 'default');
    EEG = eeg_checkset( EEG );
    
    % Let's fit one dipole to better visualize if the IC is neural or not
    EEG = pop_dipfit_settings( EEG, 'hdmfile','C:\\Program Files\\MATLAB\\R2020b\\toolbox\\eeglab2020_0\\plugins\\dipfit3.4\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','C:\\Program Files\\MATLAB\\R2020b\\toolbox\\eeglab2020_0\\plugins\\dipfit3.4\\standard_BEM\\standard_mri.mat','chanfile','C:\\Program Files\\MATLAB\\R2020b\\toolbox\\eeglab2020_0\\plugins\\dipfit3.4\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0 0 0 0 0 -1.5708 1 1 1] ,'chansel',[1:61] );
    EEG = eeg_checkset( EEG );
    EEG = pop_multifit(EEG, size(EEG.icaweights,1) ,'threshold',100,'plotopt',{'normlen','on'});
    EEG = eeg_checkset( EEG );
    
    % Save data 
    pop_saveset_with_mkdir( EEG, FILENAME_ICA, FILEPATH_OUTPUT);
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
end