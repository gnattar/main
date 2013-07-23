
function plot_dist_sessions(anm,numSess)
% anm = 'gr199201';
% numSess = 8;
 sc = get(0,'ScreenSize');
figure('position', [1000, sc(4)/10-100, sc(3)*3/10, sc(4)*3/4], 'color','w');
suptitle ([anm 'Normalized Thetaenv_distributions from Sessions']);
for i = 1:numSess
    dist = wSigSummary{1, i}.nogo_thetaenv_dist{1, 1}{1, 1};
    bins = wSigSummary{1, i}.nogo_thetaenv_bins{1, 1}{1, 1} ;
   subplot(numSess,1,i);set(gcf,'DefaultAxesColorOrder',copper(size(dist,1)));
   plot( bins,dist','linewidth',1);set(gca,'XTick',[-50 : 5 :50]);
   vline(wSigSummary{1, 2}.nogo_thetaenv_biased_barpos{1, 1}{1},'r');
   vline(wSigSummary{1, 2}.nogo_thetaenv_baseline_barpos{1, 1}{1},'k');
     
end

xlabel ('Theta env (deg)');
fnam = [ anm 'Norm_Theta_dist'];
saveas(gcf,[pwd,filesep,fnam],'tif');

figure('position', [1000, sc(4)/10-100, sc(3)*3/10, sc(4)*3/4], 'color','w');
c =copper(i);suptitle ([anm 'Cumulative Thetaenv_distributions from Sessions']);
for i = 1:numSess
    bins = wSigSummary{1, i}.nogo_thetaenv_bins{1, 1}{1, 1} ;
   temp = wSigSummary{1, i}.nogo_thetaenv_trials{1, 1}{1, 1};
	temp = reshape(temp,1,prod(size(temp)));
    dist=histnorm(temp,bins);
   subplot(numSess,1,i);
   plot( bins,dist','linewidth',2,'color',c(i,:));set(gca,'XTick',[-50 : 5 :50]);
   vline(wSigSummary{1, 2}.nogo_thetaenv_biased_barpos{1, 1}{1},'r');
    vline(wSigSummary{1, 2}.nogo_thetaenv_baseline_barpos{1, 1}{1},'k');
 
end

xlabel ('Theta env (deg)');
fnam = [ anm 'Cumm_Norm_Theta_dist'];
saveas(gcf,[pwd,filesep,fnam],'tif');
