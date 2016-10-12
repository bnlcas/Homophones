function ERPs = make_homophone_task_erps_gen(evnt) %, exp_dat)
%% This fucntion is designed to form ERPs related to the prime and 
% target stimuli given in the homophone prime experiment
%
% Inputs:
%
% evnt - 1xNUM_TRIALS structured array - contains the names and paths of
% the stimuli data
%
% exp_dat - contains record of primes and stimuli for a test block
%
% Output:
%
% ERPs - structured array contains the following fields
%   ecog_primes - ecog erps of primes
%   ecog_targets - ecog erps of targets
%   is_related_dominant - boolean, if prime is for dominant sense of target
%   is_related_subordinant - boolean, if prime is for subordinant sense
%   time_axis - time axis of erps
%   target_names - list of target names
%   prime_names - list of prime names
%   target_duration - duration (in seconds) of target
%   prime_duration - duration (in seconds) of prime

%% Sort evnt such that the sequencing makes sense:
start_times = zeros(1,length(evnt));
for i = 1:length(evnt) 
   start_times(i) = evnt(i).StartTime;
end
[~,order] = sort(start_times);
evnt = evnt(order);

%% get list of tartget/prime names/starttimes
num_trials = length(evnt)/2;
target_names = cell(num_trials,1);
target_starttimes = zeros(num_trials, 1);
prime_names = cell(num_trials, 1);
data_block = cell(num_trials,1);
prime_starttimes = zeros(num_trials, 1);
target_duration = zeros(num_trials,1);
prime_duration = zeros(num_trials,1);


for i = 1:num_trials
    prime_ind = 2*i-1;
    target_ind = 2*i;
    
    prime_names{i} = evnt(prime_ind).name;
    prime_starttimes(i) = evnt(prime_ind).StartTime;
    prime_duration(i) = evnt(prime_ind).StopTime - evnt(prime_ind).StartTime;
    
    target_names{i} = evnt(target_ind).name;
    target_starttimes(i) = evnt(target_ind).StartTime;
    target_duration(i) = evnt(target_ind).StopTime - evnt(target_ind).StartTime;
end

%% find out the association between prime and target
% load list of words
% loads homophone_list_full
% col 1 - dominant sense, col 2 - subordinant sense, col 3 - homophone
dpath = strsplit(evnt(1).dpath, filesep);
base_dir = strjoin(dpath(1:(end-2)),filesep);
%load('/Users/changlab/Documents/changrepo/matlab/analysis/Homophones/Homophone_Prime_Exp_6/WordLists/homophone_list_full.mat')
load([base_dir filesep 'behavior' filesep 'homophone_list_full.mat']);



is_related_dominant = false(num_trials,1);
is_related_subordinant = false(num_trials, 1);
for i = 1:num_trials
    is_homo = strcmpi(homophone_list_full(:,1), target_names{i});
    is_p1 = strcmpi(homophone_list_full(:,2), prime_names{i});
    is_p2 = strcmpi(homophone_list_full(:,3), prime_names{i});

    is_related_dominant(i) = (sum(is_homo & is_p1) > 0);
    is_related_subordinant(i) = (sum(is_homo & is_p2) > 0); 
end

%% Generate ECoG ERPs
%subj_dir = '/Users/changlab/Documents/data/EC123/data';
tmp = strsplit(evnt(1).dpath, filesep);
subj_dir = strjoin(tmp(1:(end-1)), filesep);
%evnt = convert_dpaths_local(evnt)

time_win = [-2 4]; % set up time range for each erp
zscore_win = [-0.5 0]; % time range prior to prime stimulus onset for z-scoring
%fs = 400;
fs = 100;
time_dim = (time_win(2)-time_win(1))*fs;
ecog_primes = [];
ecog_targets = [];
is_good_trial = true(num_trials,1);
for i = 1:num_trials
    data_block{i} = evnt(i).block;
    data_block_num = evnt(i).block;
    data_subj = evnt(i).subject;
    ecog_hg_dpath = [subj_dir filesep data_subj '_' data_block_num];
    
    zscore_range = zscore_win + [prime_starttimes(i) prime_starttimes(i)];
    prime_range = time_win + [prime_starttimes(i) prime_starttimes(i)];
    target_range = time_win + [target_starttimes(i) target_starttimes(i)];
    
    %% load data:
    ecog_z = load_ecog_data(ecog_hg_dpath, zscore_range); % for z-scoring 256x51
    ecog_p = load_ecog_data(ecog_hg_dpath, prime_range); % for prime 256x601
    ecog_t = load_ecog_data(ecog_hg_dpath, target_range); % for target 256x601
    
    % zscore ecog to before trial interval
    ecog_p.data = gdivide(gsubtract(ecog_p.data, mean(ecog_z.data,2)), std(ecog_z.data,[],2));
    ecog_t.data = gdivide(gsubtract(ecog_t.data, mean(ecog_z.data,2)), std(ecog_z.data,[],2));
    
    %% Only include good times:
    z_range_bad = contains_badTimes(ecog_z.badTimeSegments, zscore_range); % employs helper function - contains_badTimes
    p_range_bad = contains_badTimes(ecog_p.badTimeSegments, prime_range);
    t_range_bad = contains_badTimes(ecog_t.badTimeSegments, target_range); 
    is_good_trial(i) = ~z_range_bad &  ~p_range_bad & ~t_range_bad;
    
    %% add ecog to main data:
    ecog_primes = cat(3, ecog_primes, ecog_p.data(:,1:time_dim));
    ecog_targets = cat(3, ecog_targets, ecog_t.data(:,1:time_dim));
end
time_axis = linspace(time_win(1), time_win(2), time_dim);


%% Assemble ERP structure:
ERPs.ecog_primes = ecog_primes;
ERPs.ecog_targets = ecog_targets;
ERPs.time_axis = time_axis;
ERPs.is_related_dominant = is_related_dominant;
ERPs.is_related_subordinant = is_related_subordinant;
ERPs.target_names = target_names;
ERPs.prime_names = prime_names;
ERPs.BadChans = ecog_p.bad_chans;
ERPs.badTimeSegments = ecog_p.badTimeSegments;
ERPs.is_good_trial = is_good_trial;
ERPs.prime_duration = prime_duration;
ERPs.target_duration = target_duration;
ERPs.prime_starttimes = prime_starttimes;
ERPs.target_starttimes = target_starttimes;
ERPs.data_block = data_block;


end


%% Helper Functions:

%%
function ecog = load_ecog_data(dpath, time_range)
%% helper function, creates a structured array of ecog data
% input data_path, timerange (in seconds)
%
% output - ecog, fields are:
%   data - 256 x timepts matrix
%   bad_channels - list of bad channels
%   bad_timepts - list of bad_time
hg_dir = [dpath filesep 'HilbAA_70to150_8band'];
art_dir = [dpath filesep 'Artifacts'];

%% load artifacts
load([art_dir filesep 'badTimeSegments.mat']);
bad_chans = textread([art_dir filesep 'badChannels.txt']);

%% Load ecog data
blocks = 4;
ch_perblock = 64;
data = [];
for i = 1:blocks
    for j = 1:ch_perblock
        [ecog_ch, fs] = readhtk([hg_dir filesep 'Wav' num2str(i) num2str(j) '.htk'], 1000*time_range);
        % Take Mean of 8 bands if not already meaned.
        if size(ecog_ch,1) == 8 % take the average of the 8 bands for the high gamma
            ecog_ch = mean(ecog_ch,1);
        end
        if abs(fs - 400) < 0.001 % downsample if not already to 100 hz
            ecog_ch = resample(ecog_ch,1,4);
        end
        data = [data; ecog_ch];
    end
end
ecog.data = data;
ecog.badTimeSegments = badTimeSegments;
ecog.bad_chans = bad_chans;

end

%% 
function time_range_bad = contains_badTimes(bad_Segments, time_range); % employs helper function - contains_badTimes
%% Takes a nx2 matrix where each row is an interval of bad time, and 
% a 1x2 matrix demaracting a region of interest in returns true iff the ROI
% overlaps any of the 1x2 rois in the badTimeSegments matrix
time_range_bad = false;
for i = 1:size(bad_Segments,1)
    bad_segment = bad_Segments(i,:);  
% % %     w1 = abs(diff(bad_segment)); % span of the bad time segment
% % %     w2 = abs(diff(time_range)); % span of the roi
% % %     if (max(bad_segment(2),time_range(2)) - min(bad_segment(1),time_range(1))) < (w1+w2);
% % %         time_range_bad = true
% % %     end

    if min(bad_segment(2),time_range(2)) >= max(bad_segment(1), time_range(1))
        time_range_bad = true;
    end

end

end

 