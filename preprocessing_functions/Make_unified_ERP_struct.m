function ERPs = Make_unified_ERP_struct(evnt_dir)
%% This functions is designed to load each of the evnt structures for a given Patient
% and generate ERPs for each block and then combine these into a single
% block.
% This requires the existence of a set of 1xn event structures with timings
% on the stimuli for each block.
%
% Inputs:
% evnt_dir - folder containing the processed evnt files for each block of
% testing
% 
%
%current_dir = pwd;
%evnt_dir = '/Users/changlab/Documents/data/EC123/events';
%cd(evnt_dir) % this actually makes some things simple?
event_files = dir(evnt_dir); % list of all possible event structures
is_file = true(size(event_files));
for i = 1:length(event_files)
    file_name = event_files(i).name;
    if ~strcmpi(file_name(1),'B');
        is_file(i) = false;
    end
end
event_files(~is_file) = [];
%event_files(1:3) =[]; % fomatting issue.
%event_files(4:5) = [];


%% establish general ERPs structure:
ERPs.ecog_primes = [];          % 256xtimeptsxtrials prime locked erps
ERPs.ecog_targets = [];         % 256xtimeptsxtrials target locked erps
ERPs.time_axis = [];            % time axis for ecog primes and ecog targets
ERPs.is_related_dominant = [];  % 1 if target is related and dominant prime
ERPs.is_related_subordinant = [];   % 1 is target is related and suborinate prime
ERPs.target_names = [];         % cell array listing target homophones
ERPs.prime_names = [];
ERPs.BadChans = [];
ERPs.badTimeSegments = [];
ERPs.is_good_trial = [];
ERPs.prime_duration = [];
ERPs.target_duration = [];
ERPs.reaction_times = [];
ERPs.repeat_type = [];
ERPs.prime_starttimes = [];
ERPs.target_starttimes = [];
ERPs.data_block = [];


for k = 1:length(event_files)
    load([evnt_dir filesep event_files(k).name])
    % should be titled evnt by default
    
    %% sort event times:
    start_times = zeros(1,length(evnt));
    for j = 1:length(evnt)
        start_times(j) = evnt(j).StartTime;
    end
    [~,order] = sort(start_times);
    evnt = evnt(order);
    
    %% EXCEPTION HANDLING: (should be dropped in a later release - clean data before)
% % %     % drop missing pair in event 22 (would be nice to investigate this)
% % %     if strcmpi(event_files(k).name, 'B22_evnt.mat')
% % %         evnt(53) = [];
% % %     end
    %% Attempt to find the reaction times;
    reaction_times = get_response_times(evnt);
    repeat_type = get_prompt_repeat_type(evnt);
%    isi = [];
%    for i = 1:(length(evnt)/2) isi(i) = 1000*(evnt(2*i).StartTime - evnt(2*i-1).StopTime); end
%    isi_full = [isi_full isi];
    ERPs_block = make_homophone_task_erps_gen(evnt);
    
    ERPs.ecog_primes = cat(3, ERPs.ecog_primes, ERPs_block.ecog_primes);
    ERPs.ecog_targets = cat(3, ERPs.ecog_targets, ERPs_block.ecog_targets);
    
    ERPs.time_axis = ERPs_block.time_axis;
    ERPs.is_related_dominant = [ERPs.is_related_dominant; ERPs_block.is_related_dominant];
    ERPs.is_related_subordinant = [ERPs.is_related_subordinant; ERPs_block.is_related_subordinant];
    ERPs.target_names = [ERPs.target_names; ERPs_block.target_names];
    ERPs.prime_names = [ERPs.prime_names; ERPs_block.prime_names];
    ERPs.BadChans = [ERPs.BadChans, ERPs_block.BadChans];
    ERPs.badTimeSegments = [ERPs.badTimeSegments; ERPs_block.badTimeSegments];
    ERPs.is_good_trial = [ERPs.is_good_trial; ERPs_block.is_good_trial];
    ERPs.prime_duration = [ERPs.prime_duration; ERPs_block.prime_duration];
    ERPs.target_duration = [ERPs.target_duration; ERPs_block.target_duration];
    ERPs.reaction_times = [ERPs.reaction_times; reaction_times];
    ERPs.repeat_type = [ERPs.repeat_type; repeat_type];
    ERPs.prime_starttimes = [ERPs.prime_starttimes;  ERPs_block.prime_starttimes];
    ERPs.target_starttimes = [ERPs.target_starttimes; ERPs_block.target_starttimes];
    ERPs.data_block = [ERPs.data_block; ERPs_block.data_block];
end

ERPs.BadChans = unique(ERPs.BadChans);


%% Remove CH 1:64 and resected electrodes
% Exceptional situation
% ERPs.BadChans = unique([1:64 ERPs.BadChans]);
% ERPs.BadChans = unique([106:108, 122:124, 138:140,  ERPs.BadChans]);

%cd(current_dir); % return to original working directory
end