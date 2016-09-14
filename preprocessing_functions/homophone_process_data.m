function [] = homophone_process_data(subj, root_dir, blocks);
% this function finds the experimental data and audio recording for 
% the homophone experiment and generates evnt structures for the data.
%
% creates a directory of for the events created
%
% EX: 
% subj = 'EC125';
%rootdir = '/Users/changlab/Documents/data/EC125';
% blocks = {'B10'};

anin_ch = [2]; % [2 2 3 3 3 1 1 2];



write_wav_files_flag = 0;
combine_behavioral_blocks_flag = 1;
combine_ecog_blocks_flag = 0;
detect_events_flag = 1;
check_events_flag = 0;
    
%% WRITE WAV FILES

if write_wav_files_flag
    for i = 1:length(blocks)
        fprintf('Converting ANIN files for %s_B%d\n',subj,blocks{i});
        for j = 1:4
            [d,sf] = readhtk([rootdir '/dura/data_store1/human/prcsd_data/' subj ...
                '/' subj '_' blocks{i} '/Analog/ANIN' num2str(j) '.htk']);
            
            audiowrite([rootdir '/dura/data_store1/human/prcsd_data/' subj ...
                '/' subj '_' blocks{i} '/Analog/ANIN' num2str(j) '.wav'],d,round(sf));
        end
    end
end

%% COMBINE BEHAVIORAL BLOCKS

if combine_behavioral_blocks_flag
    behavior_dir = [root_dir filesep 'behavior'];
    combine_behavior_blocks(behavior_dir, false);
end

%% COMBINE ECOG BLOCKS

if combine_ecog_blocks_flag
    blocks_to_combine = {'B66','B67'};
    
    for j = 1:4
        fprintf('Writing ANIN%d\n',j);
        [d1,sf1] = readhtk([rootdir '/dura/data_store1/human/prcsd_data/' subj ...
            '/' subj '_' blocks_to_combine{1} '/Analog/ANIN' num2str(j) '.htk']);
        [d2,sf2] = readhtk([rootdir '/dura/data_store1/human/prcsd_data/' subj ...
            '/' subj '_' blocks_to_combine{2} '/Analog/ANIN' num2str(j) '.htk']);
        writehtk([rootdir '/dura/data_store1/human/prcsd_data/' subj ...
            '/' subj '_' blocks_to_combine{1} '_full/Analog/ANIN' num2str(j) '.htk'],...
            [d1 d2],sf1);
        audiowrite([rootdir '/dura/data_store1/human/prcsd_data/' subj ...
            '/' subj '_' blocks_to_combine{1} '_full/Analog/ANIN' num2str(j) '.wav'],[d1 d2],round(sf1));
    end
    
    for i = 1:4
        for j = 1:64
            fprintf('Writing Wav%d%d\n',i,j);
            [d1,sf1] = readhtk([rootdir '/dura/data_store1/human/prcsd_data/' subj ...
                '/' subj '_' blocks_to_combine{1} '/HilbAA_70to150_8band/Wav' num2str(i) num2str(j) '.htk']);
            [d2,sf2] = readhtk([rootdir '/dura/data_store1/human/prcsd_data/' subj ...
                '/' subj '_' blocks_to_combine{2} '/HilbAA_70to150_8band/Wav' num2str(i) num2str(j) '.htk']);
            writehtk([rootdir '/dura/data_store1/human/prcsd_data/' subj ...
                '/' subj '_' blocks_to_combine{1} '_full/HilbAA_70to150_8band/Wav' num2str(i) num2str(j) '.htk'],...
                [d1 d2],sf1);
            
            [d1,sf1] = readhtk([rootdir '/dura/data_store1/human/prcsd_data/' subj ...
                '/' subj '_' blocks_to_combine{1} '/RawHTK/Wav' num2str(i) num2str(j) '.htk']);
            [d2,sf2] = readhtk([rootdir '/dura/data_store1/human/prcsd_data/' subj ...
                '/' subj '_' blocks_to_combine{2} '/RawHTK/Wav' num2str(i) num2str(j) '.htk']);
            writehtk([rootdir '/dura/data_store1/human/prcsd_data/' subj ...
                '/' subj '_' blocks_to_combine{1} '_full/RawHTK/Wav' num2str(i) num2str(j) '.htk'],...
                [d1 d2],sf1);
            
        end
    end
end

%% DETECT EVENTS

if detect_events_flag
    
    %evnt = [];
    for i = 1:length(blocks)
        evnt = [];
        block_tmp = {[subj '_' blocks{i}]};
        load([root_dir, filesep 'behavior' filesep 'Homophone_Exp_' subj '_' blocks{i} '.mat']);
        audio_dat_dir = [root_dir filesep 'ecog_data'];
        stim_files = [root_dir filesep 'Event_Audio'];
        stim_names = [exp_dat.primes; exp_dat.targets];
        ev = DetectEventsQuick_phonrest(audio_dat_dir, stim_files, block_tmp, stim_names, stim_names, anin_ch, [], false);

        %load([rootdir '/data/' subj '/homophone/behavior/' subj '_' blocks{i} '.mat']);
%         ev = DetectEventsQuick_phonrest([rootdir '/dura/data_store1/human/prcsd_data/' subj],...
%             [rootdir '/tasks/homophone/findEvents'],...
%             {[subj '_' blocks{i}]},...
%             [exp_dat.primes ; exp_dat.targets],...
%             [exp_dat.primes ; exp_dat.targets],...
%             anin_ch(i),...
%             [],...
%             1);
            
        figure;
        plot([ev.confidence]);
        hold on;
        
        fprintf('Low confidence items:\n');
        fprintf('%s\t\n',strjoin({ev(find([ev.confidence] <= 0.96)).name}));
        fprintf('%2.2f\t',[ev(find([ev.confidence] <= 0.96)).confidence]);
        
        
        evnt = [evnt ev];
        event_dir = [root_dir '/events'];
        if ~exist(event_dir,'dir')
            mkdir(event_dir);
        end
        save([event_dir filesep blocks{i} '_evnt.mat'],'evnt');
    end
    
    % CHECK STIMULI VISUALLY
    if check_events_flag
        uniqueStims = unique({evnt.name});
        
        figure;
        for i = 1:length(uniqueStims)
            trls = find(strcmpi({evnt.name},uniqueStims{i}));
            nTrials(i) = length(trls);
            subplot(4,3,i);
            clear X
            for j = 1:length(trls)
                X(j,:) = x(evnt(trls(j)).StartTime*sf:evnt(trls(j)).StopTime*sf);
                %         plot(x(evnt(trls(j)).StartTime*sf:evnt(trls(j)).StopTime*sf));
                %         hold on;
            end
            r = xcorr(X','coeff');
            imagesc(reshape(max(r),sqrt(size(r,2)),sqrt(size(r,2))));
            colorbar;
            
            title(uniqueStims{i});
            %     sound(x(evnt(trls(j)).StartTime*sf:evnt(trls(j)).StopTime*sf),sf);
            %     waitforbuttonpress;
        end
    else
        fprintf('Skipping event checking....\n');
    end
else
    fprintf('Skipping event detection....\n');
end
