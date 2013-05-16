function plot_roiSignals(obj,fov,rois,roislist,tag_trialtypes,trialtypes,sfx)
% plot signals arranged by rois : to check roi selection in fovs
roisperfig = 5;
fovname = ['fov ' fov 'rois ' roislist]; 
frametime=obj.FrameTime;
rois_trials  = arrayfun(@(x) x.dff, obj,'uniformoutput',false);
numtrials = length(rois_trials);
numrois = size(rois_trials{1},1);
numframes =size(rois_trials{1},2);
ts = frametime:frametime:frametime*numframes;
temprois = zeros(numrois,numframes,numtrials);
for i= 1:numtrials
        tempmat = zeros(size(rois_trials,1),size(rois_trials,2));
        tempmat = rois_trials{i};
      if (size(tempmat,2)>size(temprois,2))
       temprois(:,:,i) = tempmat(1:numrois,1:numframes);
      else
        temprois(:,1:size(tempmat,2),i) = tempmat(:,:);
      end
      
end
cscale=[0 300];

if(tag_trialtypes ==1)
    temp = permute(temprois,[3,2,1]);
    newrois=zeros(size(temp,1),size(temp,2)+1,size(temp,3));
    newrois(:,1:size(temp,2),:) = temp;
    temp2=repmat(trialtypes,1,size(temp,3));
    temp2 = reshape(temp2,numtrials,1,numrois);
    temp2=temp2*(1/length(unique(trialtypes))) *cscale(1,2);
    temp2 = [temp2 temp2 temp2 temp2 temp2];
    newrois(:,size(temp,2)+1:size(temp,2)+5,:) = temp2;
else 
    newrois =permute(temprois,[3,2,1]);
end
numcolstoplot = 1+length(unique(trialtypes))*2;
dt = ts(length(ts))-ts(length(ts)-1);
roicount = 1;
count =1;
sc = get(0,'ScreenSize');
% h1 = figure('position', [1000, sc(4)/10-100, sc(3)*3/10, sc(4)*3/4], 'color','w');
h1 = figure('position', [1000, sc(4)/10-100, sc(3), sc(4)*1/2], 'color','w');
rois_name_tag = '';
    for i = 1:length(rois)

            rois_name_tag = [rois_name_tag,num2str(rois(i)),','];
            %plot im
            subplot(roisperfig,numcolstoplot,count);
            count=count+1;
            imagesc([ts ts(length(ts))+dt*1:dt:ts(length(ts))+dt*5],1:numtrials,newrois(:,:,rois(i)));caxis(cscale);colormap(jet);
            vline([.5 1 1.5 2 2.5 ],'k-');


            % plot traces
    % % %          if(tag_trialtypes ==1)         
                types= unique(trialtypes);          
                col = [0 0 1; 0 .5 1; 0 1 0;1 .6 0;1 0 0; .5 0 0  ];
                for k = 1:length(types)
                    trials_ktype=(find(trialtypes==types(k)));  
                    subplot(roisperfig,numcolstoplot, count);
                    count = count+1;
                    for j=1:length(trials_ktype)           
                        plot(ts ,newrois(trials_ktype(j),1:length(ts),rois(i))','color',col(types(k),:),'linewidth',1);           
                     hold on;  
                    end   
                    axis([0 ts(length(ts)) -200 800]);set(gca,'YMinorTick','on','YTick', -200:200:800);
                    vline([ .5 1 1.5 2 2.5],'k-');
                end

    % % %          else
    % % %             col=linspace(0,1,numtrials);
    % % %             for j=1:numtrials         
    % % %                  plot(ts,newrois(j,:,rois(i))','color',[.7 0 0]*col(j),'linewidth',.5);
    % % %                  axis([0 ts(length(ts)) -50 700]);set(gca,'YMinorTick','on','YTick', 0-50:200:700);
    % % %                  hold on;
    % % %             end
    % % %             hold off;  
    % % %             vline([.5 ],'r-');
    % % %             vline([1 1.5 2 2.5],'k:');
    % % %          end



            % plot trace averages
                types= unique(trialtypes);  
                temp_avg=zeros(length(types),length(ts));
                col = [0 0 1; 0 .5 1; 0 1 0;1 .6 0;1 0 0; .5 0 0  ];
                for k = 1:length(types)
                    trials_ktype=(find(trialtypes==types(k)));  
                    subplot(roisperfig,numcolstoplot, count);
                    count = count+1;
                    temp_data=newrois(trials_ktype,1:length(ts),rois(i));
                    temp_avg=mean(temp_data,1);                  
                    temp_sd=(temp_data.^2 + repmat(temp_avg.^2,size(temp_data,1),1) - 2*temp_data.*repmat(temp_avg,size(temp_data,1),1));
                    temp_sd=(mean(temp_sd,1)).^0.5;


    %               jbfill(ts,temp_avg+temp_sd,temp_avg-temp_sd,col(types(k),:),col(types(k),:),1,transparency); hold on;

                    plot(ts ,temp_avg,'color',col(types(k),:),'linewidth',1.5);           

                   axis([0 ts(length(ts)) -200 300]);set(gca,'YMinorTick','on','YTick', -200:200:300);
                    vline([.5 1 1.5 2 2.5],'k-');

                end


        if (mod(roicount,roisperfig)>0) && (roicount<length(rois))

        else
            fnam=['FOV' fov 'rois' rois_name_tag sfx '.tif'];
            set(gcf,'NextPlot','add');
            axes;
            h = title(fnam);
            set(gca,'Visible','off');
            set(h,'Visible','on'); 

             set(gcf,'PaperUnits','inches');
             set(gcf,'PaperPosition',[1 1 24 10]);
            set(gcf, 'PaperSize', [24,10]); 
            set(gcf,'PaperPositionMode','manual');

            saveas(gcf,[pwd,filesep,fnam],'tif');
            close(h1);
            if (roicount<length(rois))
    %        h1 = figure('position', [1000, sc(4)/10-100, sc(3)*3/10, sc(4)*3/4], 'color','w');
                h1 = figure('position', [1000, sc(4)/10-100, sc(3), sc(4)*1/2], 'color','w');
                count =1;
                rois_name_tag = '';
            end
        end
            roicount = roicount+1;
    end
end

function [event_detected_data] = detect_Ca_events(src_data,sampling_time)
    for i = 1:size(src_data,1)
       current_trace = src_data(i,:);
       event_filter = [zeros(round(0.1/sampling_time),1); ones(round(0.2/sampling_time),1)];
       filtersize = size(event_filter,1);
       dot_product = zeros(size(current_trace,1);
       for tau = -filtersize : 1: filtersize
            
       end
       dot_product = current_
    end
    
end

