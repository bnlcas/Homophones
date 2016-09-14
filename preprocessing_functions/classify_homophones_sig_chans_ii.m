function [class_accuracies, shuffle_accuracies, shuffle_CI_range, is_sig_anychan] = classify_homophones_sig_chans_ii(ERPs, compare_UnRel)
%% Second iteration of classifier function:
%%% The propose of this function to find significant electrodes encoding differences
% between two sets of data
%
% INPUTS: ERP data structure:
%
% VARIABLE INPUTS: (NOT YET IMPLEMENTED):
%
% 1 - Confidence Interval Width
%
% 2 - Number of Shufflings
%
% 3 - gamma setting
%
% 4 - Cross Valdation Scheme: 'KFold', 'LeaveOne Out
%
% Outputs:
% class_accuracies - (num_timepoints x num_homophones) Classification Accuracies of classifier
%
% shuffle_accuracies - shuffle accuracis
%
% shuffle_CI - confidence interval on shuffled data

%% Functions Parameters:
window_data = true;
pca_data = false;
Z_norm_timept_data = false;
save_data = false;
num_shuffles = 40;
Gamma_level = 0.6;
%gamma_levels = linspace(0,0.8,5);
if window_data
    window_size = 5;
end
if pca_data
    variance_thresh = 0.9;
end
 
if save_data
    if compare_UnRel
        dat_dir = '/Users/changlab/Documents/changrepo/matlab/analysis/Homophones/ECoG Analysis/Graphics/LDA Classification/UnrRel_826';
    else
        dat_dir = '/Users/changlab/Documents/changrepo/matlab/analysis/Homophones/ECoG Analysis/Graphics/LDA Classification/P1P2_826';
    end
end

%% Establish Priming Relations for each trial
is_good = ERPs.is_good_trial;
is_dom = (ERPs.is_related_dominant == 1) & is_good;
is_sub = (ERPs.is_related_subordinant == 1) & is_good;
is_rel = (is_dom | is_sub) & is_good;
is_unr = ~is_dom & ~is_sub & is_good;

%% Set Time Range of Analysis:
time_range = [-1 1]; % in seconds
[~,time_win(1)] = min(abs(ERPs.time_axis - time_range(1)));
[~,time_win(2)] = min(abs(ERPs.time_axis - time_range(2)));


%% Second, plot the relatedness for each homophone
homophones = unique(ERPs.target_names);
class_accuracies = zeros(length(time_win(1):time_win(2)), length(homophones));

shuffle_accuracies = zeros(length(time_win(1):time_win(2)), length(homophones));
shuffle_CI_range = zeros(length(time_win(1):time_win(2)), length(homophones));
accuracy_fig = figure;
for i = 1:length(homophones)
%for i = 1:length(gamma_levels)
%    Gamma_level = gamma_levels(i);
    homophone = homophones{i};
    is_homo = strcmpi(ERPs.target_names, homophone);

    
    %% Compare P1 P2
    if compare_UnRel % Compare Related - Unrelated Primes to Target Homophones:
        ecog_p1 = ERPs.ecog_targets(:,time_win(1):time_win(2),is_homo&is_rel);
        ecog_p2 = ERPs.ecog_targets(:,time_win(1):time_win(2), is_homo&is_unr);
    else % Compare Dominant/Subordinate Senses...
        ecog_p1 = ERPs.ecog_targets(:,time_win(1):time_win(2),is_homo&is_dom);
        ecog_p2 = ERPs.ecog_targets(:,time_win(1):time_win(2), is_homo&is_sub);
    end
    

    time_axis = ERPs.time_axis(time_win(1):time_win(2));
    
    % zero Bad Trials
    ecog_p1(ERPs.BadChans,:,:) = 0;
    ecog_p2(ERPs.BadChans,:,:) = 0;
    
    
    %% Smooth Data on Time Axis
    if window_data
        ecog_p1 = smooth_3d(ecog_p1, window_size, 2);
        ecog_p2 = smooth_3d(ecog_p2, window_size, 2);
    end

    
    
    %% Run through each time point and channel to make statistical comparison for significance
    is_sig_ttest = false(size(ecog_p1,1), size(ecog_p1,2));
    stat_tstat = zeros(size(ecog_p1,1), size(ecog_p1,2));
    
    pval = 0.01;
    %pval = 0.01/100;
    bf_correct = false; % adjust p-thresh by a boniferroni correction of # of channels
    if bf_correct
        pval = pval/(size(ecog_p1,1) - length(ERPs.BadChans));
    end
    
    for j = 1:size(ecog_p1,1)
        if (sum(j == ERPs.BadChans) == 0) % ignore badchans (t-test gives nans)
            for k = 1:size(ecog_p1,2)
                dat_p1 = squeeze(ecog_p1(j,k,:));
                dat_p2 = squeeze(ecog_p2(j,k,:));

                %% T-Test
                [is_sig_ttest(j,k),~,~,tmp] = ttest2(dat_p1, dat_p2, 'Alpha', pval);
                stat_tstat(j,k) = abs(tmp.tstat);
             end
        end
    end
    
    %% Plot the electrodes that are significant
    sig_timepts_thresh = 1; % minimum number of significant time points for a channel to be significant
    %twin_sig = [-0.3 1];% -0.5 to 1 second window about onset for significance.
    twin_sig = [0 1]; % Revision 8-29-16
    [~,t1] = min(abs(time_axis - twin_sig(1)));
    [~,t2] = min(abs(time_axis - twin_sig(2)));
    sig_window = t1:t2;
    
    sig_ch_t = find(sum(is_sig_ttest(:,sig_window ),2) >= sig_timepts_thresh);    
    % sort sig_ch_t by anatomy
    %load(['/Users/changlab/Documents/data/EC123/TDT_elecs_all.mat'])
    %anatomy = anatomy(1:256,4);
    %[~,anatomy_order]  = sort(anatomy);
    is_sig_anychan = (sum(is_sig_ttest(sig_ch_t,:),1)>0);
    
    
    %% Hand Select Electrodes?:
    %sig_ch_t = [93 220];  [147 220];
    
    p_labels = [repmat({'p1'},size(ecog_p1,3),1); repmat({'p2'},size(ecog_p2,3),1)];

    
    
    %% Generate PCA matrix for (Num_Chans x (time_pts x trials))
    if pca_data
        dat_pca = cat(3,ecog_p1(sig_ch_t,:,:), ecog_p2(sig_ch_t,:,:));
        dat_pca = reshape(dat_pca, size(dat_pca,1), size(dat_pca,2)*size(dat_pca,3)); % flatten time x trials
        dat_pca = gsubtract(dat_pca, mean(dat_pca,2)); % center data on each channel
        [pca_mat,~,~,~,exp] = pca(dat_pca');

        num_components = find(cumsum(exp)>variance_thresh*100,1);    % number of components to use in PCA projection
    else
        num_components = 0;
        pca_mat = [];
    end
    


    %% Set CV Partition:
    folds = 12;          % number of folds in K-Fold CV
    %partition = cvpartition(p_labels,'KFold',folds);
    partition = cvpartition(p_labels,'LeaveOut');
    folds = partition.NumTestSets;
    
    label_shuffles = zeros(num_shuffles, length(p_labels));
    for j = 1:num_shuffles
        label_shuffles(j,:) = randperm(length(p_labels));
    end
        
    
    %% Run Classifier 
    %% For each time point construct a classifier & Calculate accuracy
    lda_weights = zeros(length(sig_ch_t),length(time_axis));
    for j = 1:length(time_axis)
        dat1 = squeeze(ecog_p1(sig_ch_t,j,:)); % transpose necessary later
        dat2 = squeeze(ecog_p2(sig_ch_t,j,:)); % num_ch x num_trials matrix of data for the time point
        dat_comb = [dat1'; dat2'];
        
        if pca_data
            dat_comb = gsubtract(dat_comb, mean(dat_comb,1))*pca_mat(:,1:num_components);  % PCA projection of combined data   
        end
        if Z_norm_timept_data
            dat_comb = gsubtract(dat_comb, mean(dat_comb,1));
        end
        %% GENERATE CASSIFIER:     
 %       lda = fitcdiscr(dat_comb, p_labels,'Gamma',Gamma_level, 'CVPartition',partition);
 %       predicted =  kfoldPredict(lda);  % mean accuray of the k-folds for homophone i, in time point j 
 %       class_accuracies(j,i) = sum(strcmpi(predicted,p_labels))/length(p_labels);
        class_accuracies(j,i) = 0;
 
 %       lda_weighs = fitcdiscr(dat_comb, p_labels,'Gamma',Gamma_level);
 %       lda_weights(:,j) = (lda_weighs.Coeffs(1,2).Linear);
% %         for k = 1:folds
% %             train_inds =  training(partition,k);
% %             test_inds = test(partition,k);
% %             lda = fitcdiscr(dat_comb(train_inds,:), p_labels(train_inds),'Gamma',Gamma_level);
% %             accuracies(k) = sum(strcmp(predict(lda,dat_comb(test_inds,:)),p_labels(test_inds)))/length(p_labels(test_inds));
% %         end
% %         tmp = mean(accuracies)
% %         class_accuracies(j,i)

        %% Random Perm Test (For Baseline Performance
        shuffle_accuracies(j,i) = 0.5;
        CI_bounds = binon_accuracy_ci(length(p_labels),0.5,0.05);
        shuffle_CI_range(j,i) = CI_bounds(2) - CI_bounds(1);
% % % %                 accuracies = zeros(num_shuffles,1);
% % % %                 p_labels_randperm = p_labels;
% % % %                 parfor k = 1:num_shuffles
% % % %                     p_labels_randperm = p_labels(label_shuffles(k,:));
% % % %                     lda = fitcdiscr(dat_comb, p_labels_randperm,'Gamma',Gamma_level, 'CVPartition',partition);
% % % %                     predicted =  kfoldPredict(lda);  % mean accuray of the k-folds for homophone i, in time point j
% % % %                     accuracies(k) = sum(strcmpi(predicted,p_labels))/length(p_labels);
% % % %                 end
% % % %                 shuffle_accuracies(j,i) = mean(accuracies); % mean accuray of the k-folds for homophone i, in time point j
% % % %         %         n = length(p_labels)*num_shuffles;
% % % %         %         binomial_dat = accuracies*length(p_labels);
% % % %         %         binomial_cdf = cdf('Binomial',0:1:n, n, mean(accuracies));
% % % %         %         ci = 0.95;
% % % %          %       shuffle_CI_range(j,i) = (find((binomial_cdf <= (1-(1-ci)/2)),1,'last')-find((binomial_cdf>= (1-ci)/2),1))/n;
% % % %                 tmp = sort(accuracies);
% % % %                 shuffle_CI_range(j,i) = tmp(end-1)-tmp(2); % std accuracy of k-folds for ...

    end
     

     %% Plot time courses of classification accuracy (use 'LETTER' as test case)
    figure(accuracy_fig);
    subplot(3,4,i);
    %figure;
    plot(time_axis, class_accuracies(:,i),'r')
    hold on
    plot(time_axis, shuffle_accuracies(:,i), 'k--')
    plot(time_axis, shuffle_accuracies(:,i)+shuffle_CI_range(:,i)/2, 'k:')
    plot(time_axis, shuffle_accuracies(:,i)-shuffle_CI_range(:,i)/2, 'k:')
    %legend('Classifer Accuracy', 'Chance Accuracy', '95 % CI')
    ylim([0 1])
    %shadedErrorBar(time_axis, shuffle_accuracies(:,i), shuffle_CI_range(:,i)/2, 'k', 1); % mean class of k-
    shaded_patch_significant_timepoints(time_axis, is_sig_anychan, gcf)

    plot(get(gca,'XLim'), [0.5 0.5],'k--');
    plot([0 0], get(gca,'YLim'), 'k', 'LineWidth',2)
    axis tight;
    xlabel('time (s)');
    ylabel('classification accuracy');

    title({['LDA Classification Accuracy for ' upper(homophone)];['Over Time with Leave One Out CV, Gamma = ' num2str(Gamma_level)]})
    title({['LDA Accuracy for ' upper(homophone)]});
    if save_data
        savefig(gcf,[dat_dir filesep homophone '_Classifier.fig']);
    end
    
    
    plot_sig_elects = true;
    if plot_sig_elects
        plot_erp_comparison_select_electrodes(time_axis, ecog_p1, ecog_p2, sig_ch_t, is_sig_ttest)
        
        unique(ERPs.prime_names(is_homo&is_dom))
        unique(ERPs.prime_names(is_homo&is_sub))
        %ecog_unr = smooth_3d(ERPs.ecog_targets(sig_ch_t,time_win(1):time_win(2), is_homo&is_unr),window_size, 2);

% % % %         grid_dim = ceil(sqrt(length(sig_ch_t)));
% % % %         grd = reshape(1:(grid_dim^2),grid_dim,grid_dim); % grid plot...
% % % %         flatten = @(x) x(:);
% % % %         ax = gca;
% % % %         ylims = ax.YLim;
% % % %         lda_weights = max(abs(ax.YLim))*lda_weights/max(abs(lda_weights(:)));
% % % %         ylims = [min(ylims(1), min(lda_weights(:))), max(ylims(2), max(lda_weights(:)))];
% % % %       %  ylims = [min([flatten(mean(ecog_p1,3)-nansem(ecog_p1,3)); flatten(mean(ecog_p2,3)-nansem(ecog_p2,3))]), max([flatten(mean(ecog_p1,3)+nansem(ecog_p1,3)); flatten(mean(ecog_p2,3)+nansem(ecog_p1,3))])];
% % % %         for k = 1:size(lda_weights,1)
% % % %             p = plotGridPosition(grd(k), grid_dim^2, grid_dim); 
% % % %             subplot('position',p);
% % % % %            shadedErrorBar(time_axis, mean(squeeze(ecog_unr(k,:,:)),2), nansem(squeeze(ecog_unr(k,:,:)),2),[0.8 0.8 0.2],1);
% % % %              ax = gca;
% % % %              ax.YLim = ylims;
% % % %              plot(time_axis, max(abs(ax.YLim))*lda_weights(k,:)/max(abs(lda_weights(:))),'k')
% % % %              shaded_patch_significant_timepoints(time_axis, is_sig_ttest(sig_ch_t(k),:));
% % % %              plot([0 0],ylims,'k')
% % % %              ax.YLim = ylims;
% % % %         end
        if save_data
            savefig(gcf,[dat_dir filesep homophone '_ERPs.fig']);
        end
    end
    

    a = 1;
    %is_sig_any = (sum(is_sig_ttest(sig_ch_t,:),1)>0);
    %shaded_patch_significant_timepoints(time_axis, is_sig_any, gcf);
    i
    if save_data
        save([dat_dir filesep 'classifer_dat.mat'],'class_accuracies', 'shuffle_accuracies', 'shuffle_CI_range', 'is_sig_anychan');
    end
    end


% Plot average classification accuracy:
figure;
shadedErrorBar(time_axis, mean(class_accuracies,2), nansem(class_accuracies,2), 'r', 1)
%plot(repmat(time_axis,12,1)', class_accuracies)

hold on;
set(gca,'YLim', [0 1])
plot([0 0], get(gca,'YLim'), 'k', 'LineWidth',2)
plot(time_axis, mean(shuffle_accuracies,2), 'k--')
plot(time_axis, mean(shuffle_accuracies,2)+mean(shuffle_CI_range,2)/2, 'k:')
plot(time_axis, mean(shuffle_accuracies,2)-mean(shuffle_CI_range,2)/2, 'k:')
axis tight;
xlabel('time (s)');
ylabel('classification accuracy');

title({['Average LDA Classification Accuracy for All Homophones'];['Over Time with Leave One Out CV, Gamma = ' num2str(Gamma_level)]})


a = 1;



end