function evnt = DetectEventsQuick_phonrest(dtpath, wpath, expt, names, all_names, anin_to_use, anin_zero, debug_fig)
% Detect TIMIT events by matching sentence waveforms to ANIN2
%
% dtpath is path to your data
% wpath is path to the stimulus wav files,
% expt is a cell array of the block directories, 
% names is a cell array of the wav file names (without .wav) for this block only, 
% all_names is a cell array of all possible TIMIT wav files
% anin_to_use is a vector of values 1 or 2 (same length as expt) which specifies whether to
% use ANIN1 or ANIN2 for the alignment
% debug_fig =1 will plot the matches as they are calculated
%
% Usage:
%
% dtpath = '/Users/liberty/Documents/UCSF/data/timit/EC56/';
% wpath = '/Users/liberty/Documents/UCSF/matlab/TIMIT/@ECSpeech/Sounds';
% expt = {'EC56_B12','EC56_B25'};
% fnames = dir(sprintf('%s/*.wav',wpath));
% anin_to_use = [2 2]; % whether to use ANIN1 or ANIN2
% for block = 1:length(expt)
%   for iAudio = 1:length(fnames)
%      names{block}{iAudio} = strrep(fnames(iAudio).name,'.wav','');
%   end
% end
% debug_fig = 1;
% 
% evnt = DetectEventsQuick(dtpath,wpath,expt,names, all_names, anin_to_use, debug_fig )
%
% Debug plot will show the attempted "match" for each sentence - this is
% not always correct, so only matches with a correlation greater than
% corr_thresh are accepted as the appropriate event.
%
% Written May 7, 2014 by Liberty Hamilton
%
% debug_fig=1; % show debugging figure or not

new_fs = 16000;  % resampling frequency

% initialize event variables
evnt=struct();
evind = 1;

all_conf = [];

% other_fig = figure;

for blocknum = 1:length(expt) % for number of blocks
    blockstart=tic;

    if anin_to_use(blocknum) == 1
        corr_thresh = 0.2; %  or if using ANIN1 instead
    else
        corr_thresh = 0.7; % ANIN2 is cleaner, may need to change this for noisy ANIN2
    end
    
    % get data from Analog file
    dpath = [dtpath filesep expt{blocknum}];
    wname = sprintf('%s/ANIN%d.htk', dpath, anin_to_use(blocknum));
    if ~exist(wname,'file')
        wname = sprintf('%s/Analog/ANIN%d.htk', dpath, anin_to_use(blocknum));
    end
    [anin2,anin2_fs] = readhtk(wname); %load analog channel
    if ~isempty(anin_zero)
        anin2(anin_zero(1)*anin2_fs:anin_zero(2)*anin2_fs) = 0;
    end
    anin2 = resample(anin2,round(new_fs),round(anin2_fs)); %resample analog signals to 16kHz
    anin2 = anin2-mean(anin2); % subtract mean to get mean 0
    anin2 = 0.99*anin2/max(abs(anin2));
    
    % determine wav files to use for this block
    sent_wavs = names;
    
    stimnum = 1;
    while stimnum <= length(sent_wavs) % for number of audio files
        stimstart = tic;
        fprintf(1,'\nBlock %s: looking for stimulus %d out of %d\n', expt{blocknum}, stimnum, length(sent_wavs));
        
        % get sentence text
%         fid = fopen([wpath filesep sent_wavs{stimnum} '.txt']);
%         wtxt = fgetl(fid);
        tmp = regexp(sent_wavs{stimnum},'\.','split');
        wtxt = tmp{1};
%         fclose(fid);
        
        % get sentence audio
        wname = [wpath filesep sent_wavs{stimnum} '.wav'];
        [sentence, sent_fs] = audioread(wname);
        sentence = resample(sentence, new_fs, round(sent_fs));  %resample sentence audio file to 16 kHz
        if size(sentence,2) > 1
            sentence(:,2) = [];
        end
   
        fprintf('%s\n',wtxt);
        clear cc;
        % Find where this sentence audio occurs in ANIN2
        tic;
        fftlen=length(anin2);
        match_filter = cconv(flipud(sentence(:)), anin2(:), fftlen); % THIS IS MUCH FASTER THAN CONV!
        toc;
        
        % sort by maximum of the convolution, this tells you where the END
        % of the events occurred. We sort in case there is more than one
        % example of the sentence detected in the TDT signal
        [~, maxind]=max(match_filter);

        %event_starts=sort_mi(round(sort_m)==round(sort_m(end)))-length(sentence)+1; % this is where the events start
        %event_starts=sort_mi(end)-length(sentence)+1; % this is where the events start
        start_time=maxind-length(sentence)+1; % this is where the events start
        end_time=(start_time+length(sentence(:))-1);
        
        % loop through the potential start times for this sentence
        matched_sentence_segment = anin2(start_time:end_time);
        cc=corrcoef(sentence(:), matched_sentence_segment(:)); % correlation between sentence and the "match"
        all_conf = [all_conf; cc(1,2)];
        if debug_fig
            plot_match_template(sentence, matched_sentence_segment, match_filter, wtxt, cc(1,2), all_conf, corr_thresh, start_time, end_time);
        end
        if (cc(1,2)>corr_thresh)
            % save this to evnt
            evnt(evind).name = sent_wavs{stimnum}; % name of the stimulus
            evnt(evind).ind = find(strcmp(all_names, sent_wavs{stimnum})); % number of the wav file
            evnt(evind).confidence = cc(1,2); % correlation between template and signal
            evnt(evind).StartTime = start_time/new_fs; % start time in seconds
            evnt(evind).StopTime = end_time/new_fs; % stop time in seconds
            
            % zero out data from anin2
            anin2(start_time:end_time) = 0;
            
            evnt(evind).wname = wname; % name of stimulus wav file
            evnt(evind).expt = expt{blocknum}; % EC55_B6, for example
            subnametmp = strfind(expt{blocknum},'_');
            evnt(evind).subject = expt{blocknum}(1:subnametmp(1)-1); % just the subject ID
            evnt(evind).block = expt{blocknum}(regexp(expt{blocknum},'B\d+','start'):regexp(expt{blocknum},'B\d+','end')); % just the number of the block
            evnt(evind).exptind = blocknum; 
            evnt(evind).trial = stimnum; % how many events were found for this sentence
            evnt(evind).dpath = dpath; % path to the data
            fprintf(1,'Found a match for sentence (%4.3f-%4.3f), r=%3.3f\n', evnt(evind).StartTime, evnt(evind).StopTime, cc(1,2));
            evind=evind+1;
        end
        stimend=toc(stimstart);
        fprintf(1,'Time elapsed for stim: %3.3f seconds\n', stimend);
        stimnum = stimnum + 1;
    end
    timeelapsed=toc(blockstart);
    fprintf(1,'Time elapsed for block: %4.3f s\n', timeelapsed);
end

    function plot_match_template(sentence, matched_sentence_segment, match_filter, wtxt, cc, all_conf, corr_thresh, start_time, end_time)
        % plot the sentence and its best "match" from the TDT data
        if cc>corr_thresh
            clr='g';
            ttl='good';
        else
            clr='r';
            ttl='bad';
        end
        figure(1);
        subplot(4,3,[1 2 3]);
        plot(sentence,clr); axis tight;
        title(sprintf('Template: %s',wtxt));
        
        subplot(4,3,[4 5 6]);
        plot(matched_sentence_segment,clr); axis tight;
        title(sprintf('Matched TDT signal, CC=%3.3f', cc));
        
        subplot(4,3,7);
        plot(sentence, 'r'); hold on; plot(matched_sentence_segment, 'b'); hold off; axis tight;
        title(ttl); % plot sentence and match on top of each other, title is whether this is considered a good match or not
        
        subplot(4,3,8);
        plot(sentence(:), matched_sentence_segment(:),'.'); axis tight;
        xlabel('Template sentence played'); ylabel('TDT signal');
        
        subplot(4,3,9);
        edges=[-1:0.05:1];
        bar(edges,histc(all_conf(:), edges)/length(all_conf(:))); axis tight;
        ylabel('Frequency');
        xlabel('Correlation');
        title('Running correlation histogram');
        
        subplot(4,3,[10 11 12]);
        plot(match_filter); hold all;
        ff=fill([start_time start_time end_time end_time], [min(match_filter) max(match_filter) max(match_filter) min(match_filter)],'c');
        set(ff,'facealpha',0.5,'edgecolor','c','linestyle','none');
        hold off;
        title('filter'); axis tight;
        drawnow();
        %if cc>0.7
        %pause()
        %end
    end
end