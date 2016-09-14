function [] = combine_behavior_blocks(behavior_dat_dir, varargin)
% function cycles through behavioral data blocks
% and combines the subblocks for a given block of testing
% this is done by combining things that end in the same 
% 
% Input directory where homophone experiment behavior files are stored.
%
% Variable input: delete constituent files (variable true false)
if nargin > 1
    delete_constituent_files = varargin{1};
else
    delete_constituent_files = false;
end

%% Find unique blocks within the behavioral Data (Distinguished by the last _Number)
file_names = strsplit(ls(behavior_dat_dir),{'\t','\n'});
file_names(strcmpi(file_names,'')) = [];
block_tags = cell(1,length(file_names));
for i = 1:length(file_names)
    filename_split = strsplit(file_names{i},'_');
    end_tag = filename_split{end};
    if ~strcmpi(end_tag(1),'B')
        block_tags{i} = filename_split{end-1};
    else
        block_tags{i} = '';
    end   
end
unique_blocks = unique(block_tags); % Find unique blocks...
unique_blocks(strcmpi(unique_blocks,'')) = [];

% for each unique block, combine relevant events
for i = 1:length(unique_blocks)
    block_inds = find(strcmpi(block_tags, unique_blocks{i}));
    block_tag = file_names{block_inds(1)};
    tmp = strsplit(block_tag,'_'); % split by underscore
    tmp(end) = []; % kill '_subblock_num.mat'
    combined_block_name = [strjoin(tmp,'_') '.mat'];
    
    % Initailize combined block data (HARDCODE >[)
    exp_dat.subj_id = [];
    exp_dat.block_num = [];
    exp_dat.exp_start_time = [];
    exp_dat.primes = [];
    exp_dat.targets = [];
    exp_dat.cues = [];
    exp_dat.Prime_StartTimes = [];
    exp_dat.Target_StartTimes = [];
    exp_dat.Cue_StartTimes = [];
    exp_dat.Prime_Duration = [];
    exp_dat.Target_Duration = [];
    for j = 1:length(block_inds)
        tmp = load([behavior_dat_dir filesep file_names{block_inds(j)}]);
        if delete_constituent_files % Kick out the ladder you just climbed.
            delete([behavior_dat_dir filesep file_names{block_inds(j)}]);
        end
        exp_dat.subj_id = tmp.exp_dat.subj_id;
        exp_dat.block_num = tmp.exp_dat.block_num;
        exp_dat.exp_start_time = tmp.exp_dat.exp_start_time;
        exp_dat.primes = [exp_dat.primes ; tmp.exp_dat.primes];
        exp_dat.targets = [exp_dat.targets ; tmp.exp_dat.targets];
        exp_dat.cues = [exp_dat.cues ; tmp.exp_dat.cues];
        exp_dat.Prime_StartTimes = [exp_dat.Prime_StartTimes ; tmp.exp_dat.Prime_StartTimes];
        exp_dat.Target_StartTimes = [exp_dat.Target_StartTimes ; tmp.exp_dat.Target_StartTimes];
        exp_dat.Cue_StartTimes = [exp_dat.Cue_StartTimes ; tmp.exp_dat.Cue_StartTimes];
        exp_dat.Prime_Duration = [exp_dat.Prime_Duration ; tmp.exp_dat.Prime_Duration];
        exp_dat.Target_Duration = [exp_dat.Target_Duration ; tmp.exp_dat.Target_Duration];
    end

    save([behavior_dat_dir filesep combined_block_name],'exp_dat');
end

    
end
        
