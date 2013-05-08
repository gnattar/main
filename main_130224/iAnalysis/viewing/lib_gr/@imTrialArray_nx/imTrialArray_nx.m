% Class for an array of imaging trials
% NX Feb 2009
%
classdef imTrialArray_nx < handle
    
    properties
        SessionName = '';
        AnimalName = '';
        ExpDate = '';
        imTrials = {};  % take the object array of imaging trials
        % FileName = '';
        FileName_prefix = '';
        nTrials = [];
        TrialNums = [];
        nROIs = []; % total number of ROIs
        ROI_type = {};
        
        nFramesPerTrial = [];
        FrameTime =[]; % time for each frame, in ms
        nChannels = 1;
        Ca_events = {}; %
        
        
        SoloTrials = {};
        SoloTrialNums = [];
        EphusTrials = [];
        WhiskerTrials = {};
        %         SoloTrials = [];
        PPT_filename = '';
    end
    
    properties (Dependent, SetAccess = private)
        ROI_events_param
        
    end
    
    methods (Access = public)
        %%
        function obj = imTrialArray_nx(imTrialsObj, varargin)
            obj.SessionName = imTrialsObj(1).SessionName;
            obj.AnimalName = imTrialsObj(1).AnimalName;
            obj.ExpDate = imTrialsObj(1).ExpDate;
            obj.FileName_prefix = imTrialsObj(1).FileName_prefix;
            obj.nTrials = length(imTrialsObj);
            obj.TrialNums = 1:obj.nTrials;
            obj.nROIs = imTrialsObj(1).nROIs;
            obj.ROI_type = imTrialsObj(1).ROIType;
            obj.nFramesPerTrial = imTrialsObj(1).nFrames;
            obj.FrameTime = imTrialsObj(1).FrameTime;
            obj.nChannels = imTrialsObj(1).nChannel;
            for i = 1:length(imTrialsObj)
                obj.imTrials{i} = imTrialsObj(i);
                
                if ~isempty(imTrialsObj(i).behavTrial)
                    obj.SoloTrialNums(i) = imTrialsObj(i).behavTrial.trialNum;
                    obj.SoloTrials{i} = imTrialsObj(i).behavTrial;
                end
                for j = 1:obj.nROIs
                    obj.Ca_events{i,j} = imTrialsObj(i).CaTransients{j};
                end
            end
            %             obj.ROI_events_param = imTrialsObj.ROI_events_param();
            
        end
        % ************************************************************************************
        %%
        function trial_inds = trialInds_sorted_by_behav(obj)
            ind_hit = []; ind_miss=[]; ind_cr=[]; ind_fa=[];
            for i = 1:obj.nTrials
                if obj.SoloTrials{i}.trialType==1
                    if obj.SoloTrials{i}.trialCorrect==1
                        ind_hit =[ind_hit i];
                    else
                        ind_miss=[ind_miss i];
                    end
                else
                    if obj.SoloTrials{i}.trialCorrect==1
                        ind_cr=[ind_cr i];
                    else
                        ind_fa = [ind_fa i];
                    end
                end
            end
            trial_inds.go = [ind_hit ind_miss];
            trial_inds.nogo = [ind_cr ind_fa];
            trial_inds.hit = ind_hit;
            % trial_inds.hit_goPos = ind_hit_goPos;
            trial_inds.miss = ind_miss;
            trial_inds.cr = ind_cr;
            trial_inds.fa = ind_fa;
            
            % Sort go trials by object positions.
            for i = 1:length(trial_inds.go)
                goPos(i) = obj.SoloTrials{trial_inds.go(i)}.goPosition;
            end
            [goPos_sorted, inds_temp] = sort(goPos);
            trial_inds.go_pos_sort = trial_inds.go(inds_temp);
        end
        
        %%
        function ROI_overview_plot(obj, doppt, ppt_filename)
            if nargin < 2
                doppt = 0;
            end
            if nargin < 3
                ppt_filename = obj.PPT_filename;
            end
            eventsParam = obj.ROI_events_param;
            color_sc1 = [min(min(eventsParam(2).peaks(:)), min(eventsParam(3).peaks(:)))*0.8 ...
                max(max(eventsParam(2).peaks(:)), max(eventsParam(3).peaks(:)))*0.8];
            fig1 = roi_overview_color_plot(eventsParam(2).peaks,...
                'Events Peak in Stim Epoch',color_sc1); % 'Stim Epoch'
            fig2 = roi_overview_color_plot(eventsParam(3).peaks,...
                'Events Peak in Reward Epoch',color_sc1); % 'Reward Epoch'
            
            avg1 = eventsParam(2).peak_mean; se1 = eventsParam(2).peak_se;
            avg2 = eventsParam(3).peak_mean; se2 = eventsParam(3).peak_se;
            avg3 = eventsParam(1).peak_mean; se3 = eventsParam(1).peak_se;
            ystr = 'Peak dF/F (%)'; %'Event Probability'; %
            xstr = 'ROI #';
            legstr = {'Stim', 'Reward', 'Pre-Stim'};
            
            % plot trial averaged parameters for different ROIs
            fig3 = figure('Position',[900 300 500 260]); hold on;
            h(1) = errorbar(avg1, se1, 'o-');
            h(2) = errorbar(avg2, se2, 'r-o');
            h(3) = errorbar(avg3, se3, 'g-o');
            set(gca, 'FontSize', 15);
            xlim([0 length(avg1)])
            ylabel(ystr, 'FontSize', 18);
            xlabel(xstr, 'FontSize', 18);
            legend(legstr, 'FontSize', 15);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            color_sc2 = [min(min(eventsParam(2).numEvent(:)), min(eventsParam(3).numEvent(:))) ...
                max(max(eventsParam(2).numEvent(:)), max(eventsParam(3).numEvent(:)))];
            fig4 = roi_overview_color_plot(eventsParam(2).numEvent,...
                'Events Prob. in Stim Epoch',color_sc2); % 'Stim Epoch'
            fig5 = roi_overview_color_plot(eventsParam(3).numEvent,...
                'Events Prob. in Reward Epoch',color_sc2); % 'Reward Epoch'
            
            avg1 = eventsParam(2).numEvent_mean; se1 = eventsParam(2).numEvent_se;
            avg2 = eventsParam(3).numEvent_mean;  se2 = eventsParam(3).numEvent_se;
            avg3 = eventsParam(1).numEvent_mean;  se3 = eventsParam(1).numEvent_se;
            ystr = 'Event Probability'; %'Peak dF/F (%)'; %
            xstr = 'ROI #';
            legstr = {'Stim', 'Reward', 'Pre-Stim'};
            
            % plot trial averaged parameters for different ROIs
            fig6 = figure('Position',[900 300 500 260]); hold on;
            h(1) = errorbar(avg1, se1, 'o-');
            h(2) = errorbar(avg2, se2, 'r-o');
            h(3) = errorbar(avg3, se3, 'g-o');
            set(gca, 'FontSize', 15);
            xlim([0 length(avg1)])
            ylabel(ystr, 'FontSize', 18);
            xlabel(xstr, 'FontSize', 18);
            legend(legstr, 'FontSize', 15);
            if doppt == 1
                saveppt2(ppt_filename, 'figure',[fig1,fig2,fig3,fig4,fig5,fig6], 'columns',3);
            end
            
            
            function fig = roi_overview_color_plot(param, title_str, clr_sc)
                if nargin < 4
                    clr_sc = [-10 200];
                end
                cscale = [-10 200];
                fig = figure('Position',[30   240   480   580]);
                imagesc(param); colorbar; caxis(clr_sc); set(gca, 'FontSize',12);
                colormap('Hot');
                xlabel('ROI #', 'FontSize', 15);
                ylabel('Trial #', 'FontSize', 15);
                title(title_str, 'FontSize', 18);
                set(gca,'YDir','normal');
            end
            
        end
        
        %% Scatter plot, Compare Epochs for all ROIs
        function [h,p] = ROI_epoch_compare(obj,ppt_filename,trial_type)
            if nargin < 2
                ppt_filename = obj.PPT_filename;
                trial_type = 'all';
            end
            if nargin < 3
                trial_type = 'all';
            end
            inds = obj.trialInds_sorted_by_behav;
            switch trial_type
                case 'all'
                    trialInds = 1:obj.nTrials;
                case 'go'
                    trialInds = inds.go;
                case 'nogo'
                    trialInds = inds.nogo;
                case 'hit'
                    trialInds = inds.hit;
                case 'cr'
                    trialInds = inds.cr;
            end
            trialTS=(1:obj.nFramesPerTrial).*obj.FrameTime/1000;
            
            for r = 1:obj.nROIs
                %     numEvents{r} = rp_stim.numEvent(:,r);
                %     numEvents{r} = rp_rwd.numEvent(:,r);
                fig(r) = figure; hold on;
                count=0;
                peakF_stim = [];
                peakF_rwd = [];
                mdFF_stim = [];
                mdFF_rwd = [];
                for i = trialInds
                    stimOnset = obj.SoloTrials{i}.pinDescentOnsetTime;
                    stimOffset = obj.SoloTrials{i}.pinAscentOnsetTime + obj.SoloTrials{i}.waterValveDelay;
                    epochStimInd = find(trialTS>stimOnset & trialTS<stimOffset);
                    epochRewardInd = find(trialTS>stimOffset);
                    count = count+1;
                    val_epochStim = obj.imTrials(i).dff(r,epochStimInd);
                    val_epochReward = obj.imTrials(i).dff(r,epochRewardInd);
                    mdFF_stim(count) = mean(val_epochStim);
                    mdFF_rwd(count) = mean(val_epochReward);
                    peakF_stim(count) = max(val_epochStim);
                    peakF_rwd(count) = max(val_epochReward);
                end
                [h(r),p(r)] = ttest(peakF_rwd, peakF_stim);
                plot(peakF_rwd, peakF_stim, 'o');
                lim(1) = min([get(gca, 'YLim') get(gca, 'XLim')]);
                lim(2) = max([get(gca, 'YLim') get(gca, 'XLim')]);
                xlim([lim(1) lim(2)]); ylim([lim(1) lim(2)]);
                line([lim(1) lim(2)], [lim(1) lim(2)],'Color','r');
                set(gca,'FontSize', 15, 'box', 'on');
                ylabel('Peak dF/F in Stim', 'FontSize', 18);
                xlabel('Peak dF/F in Reward', 'FontSize', 18);
                title(['ROI #' num2str(r) ' [' obj.ROI_type{r} '] ' trial_type], 'FontSize', 20, 'Color', [0.5 0.2 0]);
                text(lim(1)+10, lim(2)-20, ['P=' num2str(p(r))], 'FontSize', 13, 'Color', [0 p(r)<0.05 0]);
            end
            
            saveppt2(ppt_filename,'figure',fig(1:end),'scale','halign','left','columns',ceil(sqrt(length(fig))));
        end
        
        %% batch plot color raster for each ROIs with significant modulation
        function h_figs = batch_color_raster(obj,ROIs,color_scale, doppt_flag,ppt_filename)
            if nargin < 5
                ppt_filename = obj.PPT_filename;
            end
            if nargin < 4
                doppt_flag = 0;
%                 ppt_filename = obj.PPT_filename;
            end
            trial_inds = obj.trialInds_sorted_by_behav;
            for i = 1:length(ROIs)
                figs = obj.sort_and_plot_trials_roi(ROIs(i),trial_inds,'goPos',color_scale);
                fig_tr(i) = figs(1); % plot of trials
                fig_mean(i) = figs(3); % plot of means of go and nogo trials
                delete(figs(2));
            end
            if doppt_flag == 1
                ncol = 4; % n of columns per slides in PPT file
                for i = 1: ceil(length(ROIs)/ncol)
                    ind = (i-1)*ncol+1 : min(i*ncol, length(ROIs));
                    saveppt2(ppt_filename,'figure',[fig_tr(ind) fig_mean(ind)],'scale','true',...
                        'valign','center','padding',[0 0 20 0],'columns',length(ind));
                end
            end
            h_figs = [fig_tr fig_mean];
        end
        %% Color plot of dF/F traces of trials sorted by behavior. Mark behavioral time.
        function h_figs = sort_and_plot_trials_roi(obj, ROInum, trial_inds, sorting,color_scale)
            %              Sort Ca Trial objects - Aug, 09
            tr_hit = [obj.imTrials{trial_inds.hit}];
            tr_miss = [obj.imTrials{trial_inds.miss}];
            tr_cr = [obj.imTrials{trial_inds.cr}];
            tr_fa = [obj.imTrials{trial_inds.fa}];
            tr_go = [obj.imTrials{trial_inds.go}];
            tr_nogo = [obj.imTrials{trial_inds.nogo}];
            % Sort Ca hit trials by go positions (default)
            if ~exist('sorting','var') || strcmpi(sorting,'goPos')
                goPos_hit = [];
                for i = 1:length(trial_inds.hit)
                    goPos_hit(i) = obj.SoloTrials{trial_inds.hit(i)}.goPosition;
                end;
                [goPos_hit_sort, inds_sort_by_pos] = sort(goPos_hit);
                tr_hit = tr_hit(inds_sort_by_pos);
                ind_hit = trial_inds.hit(inds_sort_by_pos);
            end
            
            %             polePos_colors = color_go_pos(obj)
            %% Label the pole position of go trial
            polePos_colors = zeros(obj.nTrials,3);
            if ~isempty(obj.SoloTrials)
                for i = 1:obj.nTrials
                    if ismember(i, trial_inds.go)
                        PolePosTrial(i) = obj.SoloTrials{i}.goPosition;
                    else
                        PolePosTrial(i) = obj.SoloTrials{i}.nogoPosition;
                    end
                end
                polePos = unique(PolePosTrial);
                nPos = numel(polePos);
                cmap = jet; cmap = cmap(15:55,:); close;
                clrs = cmap((1:nPos) * floor(size(cmap,1)/nPos),:);
                for k = 1:nPos
                    ind = find(PolePosTrial == polePos(k));
                    polePos_colors(ind,:) = repmat(clrs(k,:), numel(ind),1);
                end
            end
            
            %% Sort
            if exist('sorting','var') && strcmpi(sorting,'AnswerLick')
                AnswerLickTimes_hit = [];
                AnswerLickTimes_fa = [];
                for i = 1:length(tr_hit)
                    AnswerLickTimes_hit(i) = tr_hit(i).behavTrial.answerLickTime;
                end
                [answerT_sort_hit, inds_sort_hit] = sort(AnswerLickTimes_hit);
                tr_hit = tr_hit(inds_sort_hit);
                
                for i = 1:length(tr_fa)
                    AnswerLickTimes_fa(i) = tr_fa(i).behavTrial.answerLickTime;
                end
                [answerT_sort_fa, inds_sort_fa] = sort(AnswerLickTimes_fa);
                tr_fa = tr_fa(inds_sort_fa);
            end
            
            %

%             nTrials = length(obj);
            if ~isempty(obj.ROI_type)
                ROItype = obj.ROI_type{ROInum};
            else
                ROItype = '';
            end
            titleStr = ['ROI# ' num2str(ROInum) '(' ROItype ')' '-' obj.AnimalName '-' obj.ExpDate '-' obj.SessionName];
            color_sc = [-10 200];
            ts = (1:obj.nFramesPerTrial).*obj.FrameTime;
            if obj.FrameTime > 1
                ts = ts/1000;
            end
            if ~exist('fig1','var') || ~ishandle(fig1)
                scrsz = [1 1 1440 900]; % get(0, 'ScreenSize');
                fig1 = figure('Position', [20, 50, scrsz(3)/4+100, scrsz(4)-200], 'Color', 'w');
            else
                figure(fig1); clf;
            end;
            h_axes0 = axes('Position', [0 0 1 1], 'Visible', 'off');
            if ~isempty(tr_hit)
                h_axes(1) = axes('Position', [0.1, 0.05, 0.8, length(tr_hit)/obj.nTrials*0.85]);
                [traces_hit, hit_mean, hit_se, ts_hit] = tr_hit.get_traces_and_plot(ROInum,h_axes(1),polePos_colors);
                
                if exist('sorting','var') && strcmpi(sorting,'AnswerLick')
                    [traces_hit_align, ts1] = align_traces(CaTraces_hit, ts, answerT_sort_hit);
                    hit_mean = mean(traces_hit_align,1);
                end
                set(h_axes(1),'Box','off', 'FontSize',13,'YTickLabel','');
                
            end
            if ~isempty(tr_miss)
                h_axes(2) = axes('Position', [0.1, length(tr_hit)/obj.nTrials*0.85+0.06, 0.8,...
                    length(tr_miss)/obj.nTrials*0.85]);
                [traces_miss, miss_mean, miss_se, ts_miss] = tr_miss.get_traces_and_plot(ROInum,h_axes(2),polePos_colors);
                set(h_axes(2),'XTickLabel','','Box','off', 'FontSize',13,'YTickLabel','');
            end
            if ~isempty(tr_cr)
                h_axes(3) = axes('Position', [0.1, length([tr_hit tr_miss])/obj.nTrials*0.85+0.07, 0.8,...
                    length(tr_cr)/obj.nTrials*0.85]);
                [traces_cr, cr_mean, cr_se, ts_cr] = tr_cr.get_traces_and_plot(ROInum,h_axes(3), polePos_colors);
                set(h_axes(3),'XTickLabel','','Box','off', 'FontSize',13,'YTickLabel','');
            end
            
            if ~isempty(tr_fa)
                h_axes(4) = axes('Position', [0.1, length([tr_hit tr_miss tr_cr])/obj.nTrials*0.85+0.08,...
                    0.8, length(tr_fa)/obj.nTrials*0.85]);
                [traces_fa, fa_mean, fa_se, ts_fa] = tr_fa.get_traces_and_plot(ROInum,h_axes(4), polePos_colors);
                set(h_axes(4),'XTickLabel','','Box','off', 'FontSize',13,'YTickLabel','');
                if exist('sorting','var') && strcmpi(sorting,'AnswerLick')
                    [traces_fa_align, ts2] = align_traces(CaTraces_fa, ts, answerT_sort_fa);
                    fa_mean = mean(traces_fa_align,1);
                end
            end
            title(titleStr, 'FontSize', 18);
            if ~exist('color_scale','var') || isempty(color_scale)
                allTraces = obj.get_f_array(ROInum);
                clim(1) = round((prctile(allTraces(:),0.5))/10)*10; % round((min(allTraces))/10)*10;
                clim(2) = round((prctile(allTraces(:),99.5))/10)*10; % max(allTraces); %
            else
                clim = color_scale;
            end
            clrsc_str = ['Color Scale: [' num2str(clim(1)) ', ' num2str(clim(2)) ']'];
            axes(h_axes0); text(0.3,0.01,clrsc_str ,'FontSize',14,'Color', 'b');
            disp(clrsc_str);
            for i=1:length(h_axes)
                set(h_axes(i), 'CLim', clim);
            end
            
            
            h_figs(1) = fig1;
            
            % plot trial mean
            
            fpos = get(fig1,'Position');
            fig2 = figure('Position', [fpos(1)+100 fpos(2) fpos(3) fpos(3)/2]);
            hold on;
            if exist('sorting','var') && strcmpi(sorting,'AnswerLick')
                plot(ts1, hit_mean,'r', 'LineWidth', 1.5);
                plot(ts2, fa_mean,'k', 'LineWidth', 1.5);
                legend('Hit', 'F-A');
                set(gca,'FontSize',13,'Position',[0.14 0.24 0.805 0.68]);
                x1 = min(ts1(1),ts2(1)); x2 = max(ts1(end),ts2(end));
                xlim([x1 x2]);
                % yl = get(gca,'YLim'); ylim([-5 yl(2)]);
                set(gca,'XTick',(floor(x1):round(x2)));
                set(get(gca,'XLabel'), 'String', 'Time (sec)', 'FontSize', 18);
                set(get(gca,'YLabel'), 'String', 'mean dF/F (%)', 'FontSize',18);
            else
                % errorshade(ts, hit_mean, hit_se, 'r');
                % errorshade(ts, miss_mean, miss_se, 'b');
                % errorshade(ts, cr_mean, cr_se, 'g');
                % errorshade(ts, fa_mean, fa_se, 'm');
                if exist('hit_mean','var')
                plot(ts, hit_mean,'b', 'LineWidth', 1.5);
                end
                if exist('miss_mean','var')
                plot(ts, miss_mean, 'k', 'LineWidth', 1.5);
                end
                if exist('cr_mean','var')
                plot(ts, cr_mean,'r', 'LineWidth', 1.5);
                end
                if exist('fa_mean','var')
                plot(ts, fa_mean,'g', 'LineWidth', 1.5);
                end
                legend('Hit', 'Miss', 'C-R', 'F-A');
                set(gca,'FontSize',13,'Position',[0.14 0.24 0.805 0.68]);
                xlim([ts(1) ts(end)]);
                % yl = get(gca,'YLim'); ylim([-5 yl(2)]);
                set(gca,'XTick',(floor(ts(1)):round(ts(end))));
                set(get(gca,'XLabel'), 'String', 'Time (sec)', 'FontSize', 18);
                set(get(gca,'YLabel'), 'String', 'mean dF/F (%)', 'FontSize',18);
                
                go_mean = mean(tr_go.get_f_array(ROInum),1);
                go_se = std(tr_go.get_f_array(ROInum),0,1)./sqrt(length(tr_go));
                nogo_mean = mean(tr_nogo.get_f_array(ROInum),1);
                nogo_se = std(tr_nogo.get_f_array(ROInum),0,1)./sqrt(length(tr_nogo));
                
                fig3 = figure('Position', [fpos(1)+200 fpos(2) fpos(3) fpos(3)/2]);
                hold on;
%                 plot(ts, go_mean, 'r', 'LineWidth', 2);
%                 plot(ts, nogo_mean, 'k', 'LineWidth', 2);
                errorshade(ts, cr_mean,cr_se,[1 0 0],3, [1 .4 .3]);
                errorshade(ts, hit_mean,hit_se, [0 0 1], 3, [0.3 0.4 0.9]);
%                 legend('Go-trials', 'NoGo-trials');
                set(gca,'FontSize',13,'Position',[0.14 0.24 0.805 0.68]);
                xlim([ts(1) ts(end)]);
                % yl = get(gca,'YLim'); ylim([-5 yl(2)]);
                set(gca,'XTick',(floor(ts(1)):round(ts(end))));
                set(get(gca,'XLabel'), 'String', 'Time (sec)', 'FontSize', 18);
                set(get(gca,'YLabel'), 'String', 'mean dF/F (%)', 'FontSize',18);
                
                tr_detected = [tr_hit tr_fa];
                m_detected = mean(tr_detected.get_f_array(ROInum),1);
                se_detected = std(tr_detected.get_f_array(ROInum),0,1)/sqrt(length(tr_detected));
                tr_undetect = [tr_miss tr_cr];
                m_undetect = mean(tr_undetect.get_f_array(ROInum),1);
                se_undetect = std(tr_undetect.get_f_array(ROInum),0,1)/sqrt(length(tr_undetect));
                
                 
                fig4 = figure('Position', [fpos(1)+300 fpos(2) fpos(3) fpos(3)/2]);
                hold on;
%                 plot(ts, m_detected, 'r', 'LineWidth', 2);
%                 plot(ts, m_undetect, 'k', 'LineWidth', 2);
                errorshade(ts, m_undetect,se_undetect,[0 0 0],3, [.5 .5 .5]);
                errorshade(ts, m_detected,se_detected, [1 0 0], 3, [1 0.4 0.3]);
                legend('detected', 'undetect');
                set(gca,'FontSize',13,'Position',[0.14 0.24 0.805 0.68]);
                xlim([ts(1) ts(end)]);
                % yl = get(gca,'YLim'); ylim([-5 yl(2)]);
                set(gca,'XTick',(floor(ts(1)):round(ts(end))));
                set(get(gca,'XLabel'), 'String', 'Time (sec)', 'FontSize', 18);
                set(get(gca,'YLabel'), 'String', 'mean dF/F (%)', 'FontSize',18);
                
            end
            h_figs(2) = fig2;
            h_figs(3) = fig3;
            h_figs(4) = fig4;
            
            
            function [traces_aligned, ts_align] = align_traces(traces, ts, eventTimes)
                
                frameTime = ts(2)-ts(1);
                eventFrame = ceil(eventTimes./frameTime);
                window(1) = min(eventFrame); % ts(find(ts<min(eventTimes),1,'last'));
                window(2) = length(ts) - max(eventFrame);
                traces_aligned = [];
                for i = 1:size(traces,1)
                    inds{i} = eventFrame(i)-window(1)+1 : eventFrame(i)+window(2);
                    traces_aligned(i,:) = traces(i,inds{i});
                    ts_align = ts(inds{i})-eventTimes(i);
                end
            end
            
            
        end
        %%
        function make_2pImg_whisker_movie(obj,trialNo, wsk_mov_dir,window,save_fmt, saveFileName,FPS)
            %
            %  - NX 2009-10
            %
            imgFile2p = obj.imTrials(trialNo).FileName;
            % t_off = 5;
            frameTime = obj.FrameTime;
            if frameTime > 1
                frameTime = frameTime/1000;
            end
            
            if ~exist('saveFileName','var')
                saveFileName = obj.FileName_prefix;
            end
            
            wsk_mov_files = dir(fullfile(wsk_mov_dir,'*.seq'));
            wsk_mov_fileName = fullfile(wsk_mov_dir, wsk_mov_files(obj(trialNo).TrialNo).name);
            
            h_fig = figure('Position', [138   431   327   535]);
            ha(1) = axes('Position', [0.01 0.39 0.98 0.6]); % 320x320
            ha(2) = axes('Position', [0.01 0.01 0.98 0.38]); % 320*200
            
            im2p = imread_multi(imgFile2p,'g');
            ts2p = (1:size(im2p,3)).*frameTime;
            fr2p = find(ts2p > window(1) & ts2p <= window(2));
            count = 0;
            for i = fr2p
                axes(ha(1)); colormap(gray);
                imagesc(im2p(:,:,i),[0 400]); set(gca,'visible', 'off');
                text(3,5,['Fr# ' num2str(i)],'Color','g');
                text(400,125,[num2str(i*frameTime) ' sec'],'Color','g');
                
                t1 = (i-1)*frameTime;
                t2 = i*frameTime;
                frWsk = [round(t1/0.002) ceil(t2/0.002)];
                if frWsk(1) == 0
                    frWsk(1) = 1;
                end
                [wskImg tsWsk] = get_seq_frames(wsk_mov_fileName, frWsk, 5);
                for j = 1:size(wskImg,3)
                    axes(ha(2)); set(gca,'visible','off');
                    imshow(wskImg(:,:,j),[]);
                    text(400,30,[num2str(tsWsk(j)/1000) ' sec'],'Color','w')
                    count = count + 1;
                    F(count) = getframe(h_fig);
                    if strcmpi(save_fmt,'tif')
                        im = frame2im(F(count));
                        imwrite(im,[saveFileName '.tif'],'tif','compression','none',...
                            'writemode','append');
                    end
                end
            end
            
            if ~exist('FPS','var')
                FPS = 15;
            end
            if strcmpi(save_fmt, 'avi')
                movie2avi(F,[saveFileName '.avi'],'compression','none','fps',FPS);
            end
        end
        
        function w = get_whisker_trials(obj, WhiskerSignalTrials, whisker_ID)
            % Extract whisker data from WhiskerSignalTrials
            % INPUT: WhiskerSignalTrials, cell array containing
            %        WhiskerSignalTrial objects. If input as a filename, or
            %        [], promote user to load the data file.
            %        whisker_ID, trajectory ID starting from 0.
            % OUTPUT: struct with fields containing relevant whisker signal
            %        variables, i.e., timestamps, curvature (kappa), position (theta), touch times.
%             w = struct('ts',0,'kappa',[],'theta',[],'touch_windows',{});
            if ischar(WhiskerSignalTrials) || isempty(WhiskerSignalTrials)
                if exist(WhiskerSignalTrials,'file')
                    fn = WhiskerSignalTrials;
                else
                    fn = uigetfile('*.mat','Select Whisker Signal Array data file');
                end
                x = load(fn);
                nm = fieldnames(x);
                WhiskerSignalTrials = x.(nm{1});
            end
            wids = WhiskerSignalTrials{1}.trajectoryIDs;
            w_ind = find(whisker_ID==wids);
            for i = 1:obj.nTrials
                w.ts{i} = WhiskerSignalTrials{i}.time{w_ind};
                w.kappa{i} = WhiskerSignalTrials{i}.kappa{w_ind};
                w.theta{i} = WhiskerSignalTrials{i}.theta{w_ind};
                w.contact_events{i} = WhiskerSignalTrials{i}.contact_events{w_ind};
            end
            obj.WhiskerTrials{w_ind} = w;
        end
        
        %%
        function out = get_signal_param(obj, varargin)
            %
            % Extract signal parameters for specified ROIs, from specified
            % trials of the session.
            %
            % varargin{}, 1, ROInums
            %                 2, ROIType, 'spine', 'branch', 'trunk',  'all'.    Default, 'all'
            %                 3, trialType, 'hit', 'miss', 'cr', 'fa'.         Defualt, 'hit'
            %                 4, trialPhase, 'stim', 'reward', 'trial'.     Default, 'stim'
            %
            out = struct('ROInums', [], 'trialInds', [], 'ROItype', {'all'}, 'trialPhase','stim', 'trialType','hit',...
                'peak_dff',[], 'eventPeak',[], 'eventArea', [], 'tauDecay', [], 'fwhm', [], 'eventProb', []);
            
            if ~isempty(varargin)
                out.ROInums = varargin{1};
            end
            if length(varargin)>1 && isempty(varargin{1})
                out.ROIType = varargin{2};
            elseif isequal(length(obj.ROI_type), obj.nROIs)
                out.ROItype = obj.ROI_type(varargin{1});
            else
                out.ROItype = 'all';
            end
            if length(varargin)>2
                out.trialType = varargin{3};
            end
            if length(varargin)>3
                out.trialPhase = varargin{4};
            end
            
            ind_struct = obj.trialInds_sorted_by_behav;
            if strcmpi(out.trialType,'all')
                tr_ind = 1:obj.nTrials;
            else
                tr_ind = ind_struct.(out.trialType);
            end
            out.trialInds = tr_ind;
            
            k = find(strcmpi(out.trialPhase, {'pre_stim', 'stim', 'reward', 'trial'}));
            event_struct = obj.ROI_events_param(k);
            
            dt = obj.FrameTime/1000;
            ts = (1:obj.nFramesPerTrial).*dt;
            for i = 1:length(tr_ind)
                for j = 1:length(out.ROInums)
                    dff = obj.imTrials(tr_ind(i)).dff(out.ROInums(j),:);
                    t1 = obj.imTrials(tr_ind(i)).behavTrial.pinDescentOnsetTime;
                    t2 = obj.imTrials(tr_ind(i)).behavTrial.pinAscentOnsetTime + 0.4;
                    %                     ind_stim = find(ts>t1& ts<t2);
                    %                     ind_rwd = find(ts>t2);
                    out.peak_dff(i,j) = max(dff(ts>t1& ts<t2));
                    switch out.trialPhase
                        case 'trial'
                            out.peak_dff(i,j) = max(dff);
                        case 'stim'
                            out.peak_dff(i,j) = max(dff(ts>t1& ts<t2));
                        case 'reward'
                            out.peak_dff(i,j) = max(dff(ts>t2));
                    end
                end
            end
            out.eventPeak = event_struct.peaks(tr_ind, out.ROInums);
            out.eventArea = event_struct.areas(tr_ind, out.ROInums);
            out.tauDecay = event_struct.tauDecay(tr_ind, out.ROInums);
            out.fwhm = event_struct.fwhm(tr_ind, out.ROInums);
            out.eventProb = mean(event_struct.numEvent(tr_ind, out.ROInums));
        end
        
        function correct_lick_artifact(obj)
            
            
        end
        
        function [f_array,ts] = get_f_array(obj, roi_no, trialNums, opt)
            % Retrieve fluo time series array of particular ROIs, in the
            % form of either deltaF/F or raw fluorescence.
            if nargin < 2
                roi_no = 1: obj.nROIs;
                trialNums = 1:obj.nTrials;
                opt = 'dff'; % output deltaF/F
            end
            if nargin < 3
                trialNums = 1:obj.nTrials;
                opt = 'dff';
            end
            if nargin < 4
                opt = 'dff';
            end
            if ischar(trialNums)
                trialNums = 1:obj.nTrials;
            end
            nROIs = length(roi_no);
            f_array = zeros(length(trialNums), obj.nFramesPerTrial, nROIs);
            switch opt
                case 'raw'
                    for  i = 1:nROIs
                        for j = 1:length(trialNums)
                            f_array(j, :, i)= obj.imTrials{trialNums(j)}.f_raw(roi_no(i),:);
                        end
                    end
                case 'dff'
                    for  i = 1:nROIs
                        for j = 1:length(trialNums)
                            f_array(j, :, i)= obj.imTrials{trialNums(j)}.dff(roi_no(i),:);
                        end
                    end
            end
            ts = (1:obj.nFramesPerTrial).* obj.FrameTime/1000;
        end
        
        %%
        function out = get_mean_dff_stim_epoch(obj, roiNo, trialInd, opt)
            if nargin<4
                opt = 'mean';
            end
                
            [f,ts] = obj.get_f_array(roiNo, trialInd);
            out = [];
            for i = 1:length(trialInd)
                bp(1) = obj.SoloTrials{trialInd(i)}.pinDescentOnsetTime + 0.4;
                bp(2) = obj.SoloTrials{trialInd(i)}.pinAscentOnsetTime + 0.4;
                unitTime = obj.FrameTime; if unitTime>1, unitTime=unitTime/1000; end;
                frInd = find(ts>bp(1) & ts<bp(2));
                for j = 1:length(roiNo)
                    switch opt
                        case 'mean'
                            out(i,j) = mean(f(i, frInd, j),2);
                        case 'peak'
                            out(i,j) = max(f(i, frInd, j), [], 2);
                    end
                end
            end
        end
    end
    
    methods % Dependent property methods; cannot have attributes
        
        function ROI_events = get.ROI_events_param(obj)
            % Get the event parameters from all ROIs of all trials. And compute trial
            % ROI_events, 1x4 struct array, with each element for one of
            %                   the 4 trial phase,   'pre_stim', 'stim',  'reward', 'trial'
            %                   'pre_stim', before the pole entering
            %                   'stim', during the time the pole is presented.
            %                   'trial', get the events for the whole trial.
            %
            %
            %
            % - NX 2009
            event = struct([]);
            criteria_id = 1;
            epoch_type = {'pre_stim', 'stim', 'reward', 'trial'};
            for ii = 1:4
                ROI_events(ii).epoch = epoch_type{ii};
                pks = [];
                areas = [];
                width = nan(obj.nTrials, obj.nROIs);
                areasNorm = [];
                tauDecay = nan(obj.nTrials, obj.nROIs);
                numEvent = [];
                
                for i = 1:obj.nTrials
                    % boundary point bp
                    bp(1)=obj.SoloTrials{i} .pinDescentOnsetTime;
                    bp(2) = obj.SoloTrials{i} .pinAscentOnsetTime + obj.SoloTrials{i} .waterValveDelay;
                    unitTime = obj.FrameTime; if unitTime>1, unitTime=unitTime/1000; end;
                    switch ROI_events(ii).epoch
                        case 'pre_stim'
                            criteria_id = 1;
                            str='Pre Stim Epoch';
                            t_start=0; t_end=bp(1);
                        case 'stim'
                            str='Stim Epoch';
                            criteria_id = 2;
                            t_start=bp(1); t_end=bp(2);
                        case 'reward'
                            str='Reward Epoch';
                            criteria_id = 3;
                            t_start=bp(2); t_end=obj.nFramesPerTrial*unitTime;
                        case 'trial'
                            str='Full Trial';
                            criteria_id = 4;
                            t_start=0; t_end=obj.nFramesPerTrial*unitTime;
                    end
                    for j=1:obj.nROIs
                        event = obj.Ca_events{i,j};
                        if ~isempty(event)
                            %criteria(1:3)=[false false false];
                            for k = 1:length(event)
                                p = []; a=[]; w=[];
                                criteria(1) = event(k).time_thresh < bp(1);
                                criteria(2) = event(k).time_thresh>bp(1)&&event(k).time_thresh<bp(2); % for stim epoch
                                criteria(3) = event(k).time_thresh > bp(2); % for reward epoch
                                criteria(4) = true;
                                if criteria(criteria_id) == true
                                    p(k) = event(k).peak;
                                    a(k) = event(k).area;
                                    w(k) = event(k).fwhm;
                                    t(k) = event(k).tauDecay;
                                    [pks(i,j), ind] = max(p);
                                    areas(i,j) = a(ind);
                                    width(i,j) = w(ind);
                                    areasNorm(i,j) = max(a(ind)./w(ind));
                                    tauDecay(i,j) = t(ind); % tau of the events with largest peak
                                end
                                numEvent(i,j)= numel(p); % /(t_end-t_start);
                            end
                        else
                            pks(i,j) = NaN;
                            areas(i,j) = NaN;
                            width(i,j) = NaN;
                            tauDecay(i,j) = NaN;
                        end
                    end
                end
                
                ROI_events(ii).peaks=pks;
                ROI_events(ii).peak_mean = nanmean(pks,1);
                ROI_events(ii).peak_se = nanstd(pks,0,1)./sqrt(length(obj));
                
                ROI_events(ii).areas=areas;
                ROI_events(ii).area_mean= nanmean(areas,1);
                ROI_events(ii).area_se=nanstd(areas,0,1)/sqrt(length(obj));
                
                ROI_events(ii).fwhm=width;
                ROI_events(ii).fwhm_mean = nanmean(width,1);
                ROI_events(ii).fwhm_se= nanstd(width,0,1)/sqrt(length(obj));
                
                ROI_events(ii).areasNorm=areasNorm;
                ROI_events(ii).areasNorm_mean= nanmean(areasNorm,1);
                ROI_events(ii).areasNorm_se= nanstd(areasNorm,0,1)/sqrt(length(obj));
                
                ROI_events(ii).tauDecay = tauDecay;
                ROI_events(ii).tauDecay_mean = nanmean(tauDecay,1);
                ROI_events(ii).tauDecay_se = nanstd(tauDecay,0,1)/sqrt(length(obj));
                
                ROI_events(ii).numEvent=numEvent;
                ROI_events(ii).numEvent_mean=mean(numEvent,1);
                ROI_events(ii).numEvent_se=std(numEvent,0,1)/sqrt(length(obj));
            end
            % color plotting
        end
        
        function polePos = get_polePos(obj)
            polePos = zeros(obj.nTrials,1);
            if isempty(obj.SoloTrials)
                error('SoloTrials is empty!');
            end
            goPos = cellfun(@(x) x.goPosition, obj.SoloTrials);
            nogoPos = cellfun(@(x) x.nogoPosition, obj.SoloTrials);
            ttype = cellfun(@(x) x.trialType, obj.SoloTrials);
            polePos(ttype) = goPos(ttype);
            polePos(ttype==0) = nogoPos(ttype==0);
        end
        
        function trace_array_viewer(obj)
            imTrialViewer(obj);
        end
        
        function viewer(obj,varargin)
            %
            %   This function must be called with no arguments. Signal selection
            %       and subsequent options are then chosen through the GUI.
            %
            %   Input arguments (in varargin) are reserved for internal, recursive
            %       use of this function.
            %
            %
            %
            if nargin==1 % Called with no arguments
                objname = inputname(1); % Command-line name of this instance of a BehavTrialArray.
                h=figure('Color','white','MenuBar','none','ToolBar','Figure'); hold on;
                set(gca,'FontSize',15);
                ht = uitoolbar(h);
                a = .20:.05:0.95; b(:,:,1) = repmat(a,16,1)'; b(:,:,2) = repmat(a,16,1); b(:,:,3) = repmat(flipdim(a,2),16,1);
                bbutton = uipushtool(ht,'CData',b,'TooltipString','Back');
                fbutton = uipushtool(ht,'CData',b,'TooltipString','Forward','Separator','on');
                set(fbutton,'ClickedCallback',[objname '.viewer(''next_trial'')'])
                set(bbutton,'ClickedCallback',[objname '.viewer(''prev_trial'')'])
                hm1 = uimenu(h,'Label','ROI','Separator','on');
                      uimenu(hm1,'Label','Next','Callback',[objname '.viewer(''next_roi'')']);
                      uimenu(hm1,'Label','Previous','Callback',[objname '.viewer(''prev_roi'')']);
                hm2 = uimenu(h,'Label','Whisker ID','Separator','on');
                      uimenu(hm2,'Label','Next','Callback',[objname '.viewer(''next_wid'')']);
                      uimenu(hm2,'Label','Previous','Callback',[objname '.viewer(''prev_wid'')']);
                      
                hm3 = uimenu(h,'Label','What to Plot', 'Separator','on');
                      uimenu(hm3,'Label','dff','Checked','on','Callback',[objname,'.viewer(''plot_dff'')']);
                      uimenu(hm3,'Label','w.kappa','Checked','on','Callback',[objname,'.viewer(''plot_kappa'')']);
                      uimenu(hm3,'Label','w.theta','Checked','on','Callback',[objname,'.viewer(''plot_theta'')']);
                      uimenu(hm3,'Label','w.velocity','Checked','off','Callback',[objname,'.viewer(''plot_vel'')']);
                      
                      
                uimenu(h,'Label','Jump to trial','Separator','on','Callback',[objname '.viewer(''jumpToTrial'')']);
                
                g = struct('trialNo',1,'trialList','','dff',[],'ts',[],'roiNo',1,'wid',0,...
                    'plot_dff',1,'plot_kappa',1,'plot_theta',1,'plot_touch',1,'plot_vel',0,'w',[],'nTrials',obj.nTrials); 
                [g.dff, g.ts] = obj.get_f_array(g.roiNo);
                if isempty(obj.WhiskerTrials)
                    g.w = obj.get_whisker_trials([],g.wid);
                else
                    g.w = obj.WhiskerTrials{g.wid+1};
                end
                
                [dff ts] = obj.get_f_array(1);

            else
                g = get(gcf,'UserData');
                if isempty(g) 
                    g = struct('trialNo',1,'trialList','');
                end
                for j = 1:length(varargin);
                    argString = varargin{j};
                    switch argString
                        case 'next_trial'
                            if g.trialNo < g.nTrials
                                g.trialNo = g.trialNo + 1;
                            end
                        case 'prev_trial'
                            if g.trialNo > 1
                                g.trialNo = g.trialNo - 1;
                            end
                        case 'jumpToTrial'
                            if isempty(g.trialList)
                                nTr = obj.nTrials;
                                g.trialList = cell(1,nTr);
                                for k=1:nTr
                                    g.trialList{k} = [int2str(k) ': trialNum=' int2str(obj.TrialNums(k))];
                                end
                            end
                            [selection,ok]=listdlg('PromptString','Select a trial:','ListString',...
                                g.trialList,'SelectionMode','single');
                            if ~isempty(selection) && ok==1
                                g.trialNo = selection;
                            end
                        case 'next_roi'
                            g.roiNo = g.roiNo + 1;
                            g.dff = obj.get_f_array(g.roiNo);
                        case 'prev_roi'
                            
                            g.roiNo = g.roiNo - 1;
                            g.dff = obj.get_f_array(g.roiNo);
                        case 'plot_dff'
                            if strcmp(get(gcbo,'Checked'),'on')
                                set(gcbo,'Checked','off')
                                g.plot_dff = 0;
                            else
                                set(gcbo,'Checked','on')
                                g.plot_dff = 1;
                            end
                        case 'plot_kappa'
                            if strcmp(get(gcbo,'Checked'),'on')
                                set(gcbo,'Checked','off')
                                g.plot_kappa = 0;
                            else
                                set(gcbo,'Checked','on')
                                g.plot_kappa = 1;
                            end
                        case 'plot_theta'
                            if strcmp(get(gcbo,'Checked'),'on')
                                set(gcbo,'Checked','off')
                                g.plot_theta = 0;
                            else
                                set(gcbo,'Checked','on')
                                g.plot_theta = 1;
                            end
                        case 'plot_vel'
                            if strcmp(get(gcbo,'Checked'),'on')
                                set(gcbo,'Checked','off')
                                g.plot_vel = 0;
                            else
                                set(gcbo,'Checked','on')
                                g.plot_vel = 1;
                            end
                        case 'next_wid'
                            g.wid = g.wid + 1;
                            g.w = obj.WhiskerTrials{g.wid};
                        case 'prev_wid'
                            g.wid = g.wid - 1;
                            g.w = obj.WhiskerTrials{g.wid};
                            
                        otherwise
                            error('Invalid string argument.')
                    end
                end                
            end
            
            cla
            if g.plot_dff == 1
                g.hf = plot(g.ts,g.dff(g.trialNo,:)/100,'b-','LineWidth',2);
            end
            if g.plot_kappa == 1
                ka = g.w.kappa{g.trialNo};
                g.hka = plot(g.w.ts{g.trialNo}, medfilt1(ka)/mean(abs(ka)),'Color','k');
            end
            if g.plot_theta == 1
                th = g.w.theta{g.trialNo};
                g.hth = plot(g.w.ts{g.trialNo}, medfilt1(th)/mean(abs(th)),'r--');
            end
            if g.plot_touch == 1
                touch_events = g.w.contact_events{g.trialNo};
                if ~isempty(touch_events)
                    wsk_ts = g.w.ts{g.trialNo}; 
                    y0 = get(gca,'YLim');
                    x1 = wsk_ts(cellfun(@(x) x.inds(1), touch_events));
                    x2 = wsk_ts(cellfun(@(x) x.inds(end), touch_events));
                    x = [x1; x1; x2; x2];
                    y = repmat([y0(1); y0(2); y0(2); y0(1)],[1,size(x,2)]);
                    h_touch = patch(x,y,'k','FaceAlpha',0.2,'EdgeColor','none');
                end
            end
            if g.plot_vel == 1
                vel = diff(g.w.theta{g.trialNo})./diff(g.w.ts{g.trialNo});
                g.hv = plot(g.w.ts{g.trialNo}(2:end), medfilt1(vel,5)/mean(abs(vel)),'g-');
            end
                    
           title(['Trial ' int2str(g.trialNo) '/' int2str(obj.nTrials)]);% ', ' titleString]);
                       
            set(gcf,'UserData',g);
        end
            
    end
end
