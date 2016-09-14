function ERPs_EC123 = Make_unified_EC123_erps()
%% This functions serves as a script to load the relevant data 
current_dir = pwd;
evnt_dir = '/Users/changlab/Documents/data/EC123/events';
cd(evnt_dir)
event_files = dir;
event_files(1:2) =[];
%event_files(4:5) = [];


%% establish general ERPs structure:
ERPs_EC123.ecog_primes = [];
ERPs_EC123.ecog_targets = [];
ERPs_EC123.time_axis = [];
ERPs_EC123.is_related_dominant = [];
ERPs_EC123.is_related_subordinant = [];
ERPs_EC123.target_names = [];
ERPs_EC123.prime_names = [];
ERPs_EC123.BadChans = [];
ERPs_EC123.badTimeSegments = [];
ERPs_EC123.is_good_trial = [];
ERPs_EC123.prime_duration = [];
ERPs_EC123.target_duration = [];
ERPs_EC123.reaction_times = [];
ERPs_EC123.repeat_type = [];
ERPs_EC123.prime_starttimes = [];
ERPs_EC123.target_starttimes = [];
ERPs_EC123.data_block = [];

% isi_full = [];
for i = 1:length(event_files)
    load(event_files(i).name)
    % should be title evnt
    
    %% sort event times:
    start_times = zeros(1,length(evnt));
    for j = 1:length(evnt)
        start_times(j) = evnt(j).StartTime;
    end
    [~,order] = sort(start_times);
    evnt = evnt(order);
    
    
    % drop missing pair in event 22 (would be nice to investigate this)
    if strcmpi(event_files(i).name, 'B22_evnt.mat')
        evnt(53) = [];
    end
    %% Attempt to find the reaction times;
    reaction_times = get_response_times(evnt);
    repeat_type = get_prompt_repeat_type(evnt);
%    isi = [];
%    for i = 1:(length(evnt)/2) isi(i) = 1000*(evnt(2*i).StartTime - evnt(2*i-1).StopTime); end
%    isi_full = [isi_full isi];
    ERPs_block = make_homophone_task_erps(evnt);
    
    ERPs_EC123.ecog_primes = cat(3, ERPs_EC123.ecog_primes, ERPs_block.ecog_primes);
    ERPs_EC123.ecog_targets = cat(3, ERPs_EC123.ecog_targets, ERPs_block.ecog_targets);
    
    ERPs_EC123.time_axis = ERPs_block.time_axis;
    ERPs_EC123.is_related_dominant = [ERPs_EC123.is_related_dominant; ERPs_block.is_related_dominant];
    ERPs_EC123.is_related_subordinant = [ERPs_EC123.is_related_subordinant; ERPs_block.is_related_subordinant];
    ERPs_EC123.target_names = [ERPs_EC123.target_names; ERPs_block.target_names];
    ERPs_EC123.prime_names = [ERPs_EC123.prime_names; ERPs_block.prime_names];
    ERPs_EC123.BadChans = [ERPs_EC123.BadChans, ERPs_block.BadChans];
    ERPs_EC123.badTimeSegments = [ERPs_EC123.badTimeSegments; ERPs_block.badTimeSegments];
    ERPs_EC123.is_good_trial = [ERPs_EC123.is_good_trial; ERPs_block.is_good_trial];
    ERPs_EC123.prime_duration = [ERPs_EC123.prime_duration; ERPs_block.prime_duration];
    ERPs_EC123.target_duration = [ERPs_EC123.target_duration; ERPs_block.target_duration];
    ERPs_EC123.reaction_times = [ERPs_EC123.reaction_times; reaction_times];
    ERPs_EC123.repeat_type = [ERPs_EC123.repeat_type; repeat_type];
    ERPs_EC123.prime_starttimes = [ERPs_EC123.prime_starttimes;  ERPs_block.prime_starttimes];
    ERPs_EC123.target_starttimes = [ERPs_EC123.target_starttimes; ERPs_block.target_starttimes];
    ERPs_EC123.data_block = [ERPs_EC123.data_block; ERPs_block.data_block];
end

ERPs_EC123.BadChans = unique(ERPs_EC123.BadChans);


%% Remove CH 1:64 and resected electrodes
ERPs_EC123.BadChans = unique([1:64 ERPs_EC123.BadChans]);
ERPs_EC123.BadChans = unique([106:108, 122:124, 138:140,  ERPs_EC123.BadChans]);

cd(current_dir); % return to original working directory
end