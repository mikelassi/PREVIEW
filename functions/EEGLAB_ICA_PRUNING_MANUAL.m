% Written by: Michael Lassi
% michael.lassi@santannapisa.it

EEG = pop_loadset(PROCESSING_LIST.FILENAME_ICA{i}, PROCESSING_LIST.FILEPATH_OUTPUT{i});
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

% Threshold for automatic accepting/rejecting a component
% Don't set it less than 1/7
THRESHOLD = 0.75; % > 1 makes it completely manual

% Accept brain components over confidence threshold, reject all others
% components over threshold ("other" components are always kept to be
% revised)
ic_labels = EEG.etc.ic_classification.ICLabel;
[~, class_selection] = max(ic_labels.classifications,[],2);
manual_selection = find(all(not(ic_labels.classifications >= THRESHOLD), 2));
automatic_selection = find(any(ic_labels.classifications >= THRESHOLD, 2));
automatic_label= {ic_labels.classes{class_selection(automatic_selection)}};

% Keep in the data only brain and other when there is high confidence of
% ICLabel
to_remove = ismember(automatic_label, {ic_labels.classes{2:6}});
to_revise = automatic_selection(ismember(automatic_label, {ic_labels.classes{7}}));
automatic_removed = automatic_selection(to_remove);

% Opens window for selecting components
% To have the best visualization (with ICLabel and Dipfit enable, you have
% to modify 
EEG = pop_selectcomps(EEG,unique([manual_selection; to_revise])');

waitfor( findobj('parent', gcf, 'string', 'OK'), 'userdata');
manual_removed = find(EEG.reject.gcompreject);

% Fill additional info field to keep track of preprocessing results!
EEG.additional_info.retained_ic_var = 1 - sum(var(EEG.icaact([manual_removed'; automatic_removed],:),0,2))./ sum(var(EEG.icaact(:,:),0,2));
EEG.additional_info.removed_ic_label = groupcounts(class_selection');
EEG.additional_info.percent_ic_removed = length([manual_removed'; automatic_removed])./size(EEG.icaact,1);

% Actual removal of the components 
EEG = pop_subcomp( EEG, unique([manual_removed'; automatic_removed])', 0);

EEG = eeg_checkset( EEG );
% Final additional infos
EEG.additional_info.n_bad_channels =  length(EEG.etc.noiseDetection.reference.badChannels.all);
EEG.additional_info.bad_channels_idx = EEG.etc.noiseDetection.reference.badChannels.all;

pop_saveset_with_mkdir( EEG, FILENAME_PREPROCESSED, FILEPATH_OUTPUT);
pop_eegplot(EEG,1, 0, 0);
uiwait
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
