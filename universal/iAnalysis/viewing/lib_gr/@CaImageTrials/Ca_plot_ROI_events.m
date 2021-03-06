function ROI_events = Ca_plot_ROI_events(obj,epoch, param_plot)

% Get the event parameters from all ROIs of all trials. And comput trial
% averaging for all ROIs
% param_plot: specify which parameter to plot, peak, area, count.

% - NX 2009


event = struct([]); 
criteriaID = 1;

for i = 1:length(obj)
    % boundary point bp
    bp(1)=obj(i).behavTrial.pinDescentOnsetTime;
    bp(2) = obj(i).behavTrial.pinAscentOnsetTime + obj(i).behavTrial.waterValveDelay;
    unitTime = obj(i).FrameTime; if unitTime>1, unitTime=unitTime/1000; end;
    switch epoch
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
            t_start=bp(2); t_end=obj(i).nFrames*unitTime;
        case 'trial'
            str='Full Trial';
            criteria_id = 4;
            t_start=0; t_end=obj(i).nFrames*unitTime;
    end
    for j=1:obj(1).nROIs
        event = obj(i).CaTransients{j};
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
                    [peaks(i,j), ind] = max(p);
                    areas(i,j) = a(ind); 
                    width(i,j) = w(ind);
                    areasNorm(i,j) = max(a(ind)./w(ind)); 
                    tauDecay(i,j) = t(ind); % tau of the events with largest peak
                end
                numEvent(i,j)= numel(p); % /(t_end-t_start);
            end
        end
    end
end
% peaks(peaks==0)=NaN;
% areas(areas==0)=NaN;
% width(width==0)=NaN;
% areasNorm(areasNorm==0)=NaN;
% tauDecay(tauDecay==0)=NaN;
% for i=1:size(peaks,2)
%     NumResponse(i) = length(find(~isnan(peaks(:,i))));
% end

ROI_events.peaks=peaks;
ROI_events.peak_mean=mean(peaks,1);
ROI_events.peak_se=std(peaks,0,1)./sqrt(length(obj));
% ROI_events.peak_mean=nanmean(peaks,1);
% ROI_events.peak_se=nanstd(peaks,0,1)./sqrt(NumResponse);

ROI_events.areas=areas;
ROI_events.area_mean=mean(areas,1);
ROI_events.area_se=std(areas,0,1)/sqrt(length(obj));
% ROI_events.area_mean=nanmean(areas,1);
% ROI_events.area_se=nanstd(areas,0,1)/sqrt(NumResponse);

ROI_events.fwhm=width;
ROI_events.fwhm_mean=mean(width,1);
ROI_events.fwhm_se=std(width,0,1)/sqrt(length(obj));
% ROI_events.fwhm_mean=nanmean(width,1);
% ROI_events.fwhm_se=nanstd(width,0,1)/sqrt(NumResponse);

ROI_events.areasNorm=areasNorm;
ROI_events.areasNorm_mean=mean(areasNorm,1);
ROI_events.areasNorm_se=std(areasNorm,0,1)/sqrt(length(obj));
% ROI_events.areasNorm_mean=nanmean(areasNorm,1);
% ROI_events.areasNorm_se=nanstd(areasNorm,0,1)/sqrt(NumResponse);

ROI_events.tauDecay=tauDecay;
ROI_events.tauDecay_mean=mean(tauDecay,1);
ROI_events.tauDecay_se=std(tauDecay,0,1)/sqrt(length(obj));
% ROI_events.tauDecay_mean=nanmean(tauDecay,1);
% ROI_events.tauDecay_se=nanstd(tauDecay,0,1)/sqrt(NumResponse);

ROI_events.numEvent=numEvent;
ROI_events.numEvent_mean=mean(numEvent,1);
ROI_events.numEvent_se=std(numEvent,0,1)/sqrt(length(obj));

% color plotting
switch param_plot
    case 'peak'
        cscale = [-10 200];
        fig = figure('Position',[30   240   480   580]);
        imagesc(peaks); colorbar; caxis(cscale); set(gca, 'FontSize',12);
        %  colormap('Hot');
        xlabel('ROI #', 'FontSize', 15); ylabel('Trial #', 'FontSize', 15);
        title(['Peak Amp. of Events in ' str], 'FontSize', 18);
        set(gca,'YDir','normal');
    case 'area'
        fig = figure('Position',[50   240   480   580]);
        imagesc(areas); colorbar; set(gca, 'FontSize',12); caxis([0 250]);
        % colormap('Hot');
        xlabel('ROI #', 'FontSize', 15); ylabel('Trial #', 'FontSize', 15);
        title(['Area of Events in ' str], 'FontSize', 18);
        set(gca,'YDir','normal');
    case 'count'
        fig = figure('Position',[70   240   480   580]);
        imagesc(numEvent); colorbar; set(gca, 'FontSize',12);
        % colormap('Hot');
        xlabel('ROI #', 'FontSize', 15); ylabel('Trial #', 'FontSize', 15);
        title(['Number of Events in ' str], 'FontSize', 18);
        set(gca,'YDir','normal');
    otherwise
        return
end