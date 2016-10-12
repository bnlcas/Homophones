function is_sig_chan = find_sig_chans_homophones(ERPs, is_class1, is_class2)
%% function to determine the significant channels for subsequent analyses of the homophone task
% 
% This function will use two boolean arrays to select trials that are
% relevant for comparison and then find channels which show significant
% differnces in the TARGET LOCKED ERPs of these two sets of trials
%
% INPUTS: ERP data structure
%
% is_class1 - n_trials long boolean array that is ON for trails that are in class1
%
% is_class2 - n_trials long boolean array that is ON for class2
%
% Outputs:
% is_sig_chan - n_chanenl (256) boolean array that is true if a channel
% meet the sigificance criteria of this function
%
% Note that this only finds differences between target locked ERPs (for the
% analysis of Prime locked data, modification will be necessary).
% 
% Example: Find significant channels for comparison of related/vs unrelated
% primes.
% is_sig_chan_rel_url = find_sig_chans_homophones(ERPs, is_related, is_unrelated)
%
% Created 9-14-16 by Ben Lucas


%% Window flag - was used in most comparisions (may not be advantageous)
window_data = true;
if window_data
    window_size = 5; % smooth 5 point wide box filter
end



%% Set Time Range of Analysis:
time_range = [0 1]; % in seconds
[~,time_win(1)] = min(abs(ERPs.time_axis - time_range(1)));
[~,time_win(2)] = min(abs(ERPs.time_axis - time_range(2)));


    
%% Get ECoG for relevant Trials
ecog_c1 = ERPs.ecog_targets(:,time_win(1):time_win(2),is_class1);
ecog_c2 = ERPs.ecog_targets(:,time_win(1):time_win(2), is_class1);
    
time_axis = ERPs.time_axis(time_win(1):time_win(2));
    
% zero Bad Trials
ecog_c1(ERPs.BadChans,:,:) = 0;
ecog_c2(ERPs.BadChans,:,:) = 0;

% Smooth Data on Time Axis
if window_data
    ecog_c1 = smooth_3d(ecog_c1, window_size, 2);
    ecog_c2 = smooth_3d(ecog_c2, window_size, 2);
end


%% Run through each time point and channel to make statistical comparison for significance
is_sig_ttest = false(size(ecog_c1,1), size(ecog_c1,2));
stat_tstat = zeros(size(ecog_c1,1), size(ecog_c1,2));
    
p_thresh = 0.01;
bf_correct = false; % adjust p-thresh by a boniferroni correction of # of channels
if bf_correct
    p_thresh = p_thresh/(size(ecog_c1,1) - length(ERPs.BadChans));
end

is_sig_ttest = ttest2(ecog_c1, ecog_c2, 'Alpha', p_thresh, 'dim',3);
% % % for j = 1:size(ecog_c1,1)
% % %     if (sum(j == ERPs.BadChans) == 0) % ignore badchans (t-test gives nans)
% % %         for k = 1:size(ecog_c1,2)
% % %             dat_p1 = squeeze(ecog_c1(j,k,:));
% % %             dat_p2 = squeeze(ecog_c2(j,k,:));
% % % 
% % %             %% T-Test
% % %             test = ttest2(dat_p1, dat_p2, 'Alpha', p_thresh);
% % %             [is_sig_ttest(j,k),~,~,tmp] = ttest2(dat_p1, dat_p2, 'Alpha', p_thresh);
% % %             stat_tstat(j,k) = abs(tmp.tstat);
% % %          end
% % %     end
% % % end

%% Find Channels that show significant differences:
sig_timepts_thresh = 1; % minimum number of significant time points for a channel to be significant

is_sig_chan = (sum(is_sig_ttest,2) >= sig_timepts_thresh);    


end



%% Functions
function [out] = smooth_3d(input, smoothing, smooth_dimension);
%% Takes an 3d matrix and smooths the data in one dimension
% inputs:
% input: n x m x p matrix
% 
% smoothing: degree of smoothing
% 
% smooth_dimension: the axis to smooth the data on (1, 2, 3)
out = zeros(size(input));

if smooth_dimension == 1
    for i = 1:size(input,2)
        for j = 1:size(input,3)
            out(:,i,j) = smooth(squeeze(input(:,i,j)), smoothing);
        end
    end
elseif smooth_dimension == 2
    for i = 1:size(input,1)
        for j = 1:size(input,3)
            out(i,:,j) = smooth(squeeze(input(i,:,j)), smoothing);
        end
    end
elseif smooth_dimension == 3
    for i = 1:size(input,1)
        for j = 1:size(input,2)
            out(i,j,:) = smooth(squeeze(input(i,j,:)), smoothing);
        end
    end
end

end