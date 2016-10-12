%% Wrapper to generate ERP structure for data.

%% Example Data:
% Must list the director for the subject's data
% This folder should contain the following folders:
%
% Event_Audio - audio files for homophones
%
% behavior - contains the experiment data files from the task
% This should also contain the file 'homophone_list_full.mat'
% (This is important for determining dominant/subordinance)
%
% ecog_data - contains subfolders for each block containing the the Hilbert
% Transfored data, bad trials, bad channels and audio
% % Each block folder should be formatted
% ECXXX_BYY (Ex: EC125_B100)
% 
% % root_dir = '/Users/changlab/Documents/data/EC125';  % Subject directory described above
% % subj = 'EC125';
% % blocks = {'B46', 'B10', 'B35', 'B38', 'B40', 'B29', 'B31'};
root_dir = '/Users/changlab/Documents/data/EC131';  % Subject directory described above
subj = 'EC131';
blocks = {'B19', 'B21', 'B26', 'B31', 'B35', 'B42', 'B49'};



%% Conbine behavior data and generate homophone evnt files
homophone_process_data(subj, root_dir, blocks);


    
%% Manually Identify Mis-Identifications
% % % Fix bug on EC123_B22 (trail 53 dropped)
% % load([root_dir filesep 'events' filesep 'B31_evnt.mat']);
% % for i = 1:length(evnt) starts(i) = evnt(i).StartTime; end
% % [~,order] = sort(starts);
% % evnt = evnt(order);
% % evnt(53) = [];
% %
% % % Fix bug on EC125_B31 (trial 121 dropped)
% % load([root_dir filesep 'events' filesep 'B31_evnt.mat']);
% % for i = 1:length(evnt) starts(i) = evnt(i).StartTime; end
% % [~,order] = sort(starts);
% % evnt = evnt(order);
% % evnt(121) = [];

% % % Fix bug on EC131_B21 (trial 135 dropped)
% % load([root_dir filesep 'events' filesep 'B21_evnt.mat']);
% % for i = 1:length(evnt) starts(i) = evnt(i).StartTime; end
% % [~,order] = sort(starts);
% % evnt = evnt(order);
% % evnt(135) = [];

%% Save modified evnt files with mis-ID's removed
% % save([root_dir filesep 'events' filesep 'B31_evnt.mat']);
% % save([root_dir filesep 'events' filesep 'B21_evnt.mat']);

%% Create ERP structure from the individual blocks
evnt_dir = [root_dir filesep 'events'];
localize_evnt_files(root_dir, evnt_dir); % necessary for preprocessing
ERPs = Make_unified_ERP_struct(evnt_dir);


%% Plot Test Data
is_good = (ERPs.is_good_trial == 1);    % Trials with no bad time segments
is_sub = (ERPs.is_related_subordinant == 1) & is_good;  % subordinate sense primed
is_dom = (ERPs.is_related_dominant == 1) & is_good; % dominante sense primed
is_rel = is_sub | is_dom;                       % Related prime
is_unr = is_good & ~is_sub &~is_dom;            % Rando prime
is_mouse = is_good & strcmpi(ERPs.target_names,'mouse'); % target homophone is mouse

PlotECogGrid_std(ERPs, true, [-1 2], ERPs.ecog_targets(:,:,is_rel & is_mouse), ERPs.ecog_targets(:,:,is_unr & is_mouse))

%% Get list of significant channels for some comparisons
is_sig_chan_rel = find_sig_chans_homophones(ERPs, is_rel, is_unr);
is_sig_chan_rel_mouse = find_sig_chans_homophones(ERPs, is_rel & is_mouse, is_unr & is_mouse);



%% Save?
save_erp_struct = false;
if save_erp_struct
    save([rootdir filsep subj '_ERP_structure.mat', 'ERPs','-v7.3'); % Is freaking huge...
end


