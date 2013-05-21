function plot_roiSignals(obj,fov,rois,roislist,tag_trialtypes,trialtypes,sfx)
% plot signals arranged by rois : to check roi selection in fovs
roisperfig = 5;
fovname = ['fov ' fov 'rois ' roislist]; 
frametime=obj.FrameTime;
rois_trials  = arrayfun(@(x) x.dff, obj,'uniformoutput',false);
dKappa = cell2mat(arrayfun(@(x) x.deltaKappa{1}, obj,'uniformoutput',false)');
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
numcolstoplot = 1+length(unique(trialtypes))*2+1;
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
                    axis([0 ts(length(ts)) -200 600]);set(gca,'YMinorTick','on','YTick', -200:200:600);
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
                    [all_data,detected] = detect_Ca_events(temp_data,frametime,80);
                    detected_data= all_data(find(detected),:);
                    detected_avg=sum(detected_data ,1)./max(sum(detected,1),1);                  
                    detected_sd=(detected_data.^2 + repmat(detected_avg.^2,size(detected_data,1),1) - 2*detected_data.*repmat(detected_avg,size(detected_data,1),1));
                    detected_sd=(mean(detected_sd,1)).^0.5;


    %               jbfill(ts,temp_avg+temp_sd,temp_avg-temp_sd,col(types(k),:),col(types(k),:),1,transparency); hold on;
                    plot([frametime:frametime:length(detected_avg)*frametime] ,detected_avg,'color',col(types(k),:),'linewidth',1.5);           

                   axis([0 ts(length(ts)) -200 600]);set(gca,'YMinorTick','on','YTick', -200:200:600);

                    vline([.5 1 1.5 2 2.5],'k-');

                end
                
                
            % plot max(dFF) vs. dKappa
                
                types= unique(trialtypes);  
                    subplot(roisperfig,numcolstoplot, count);
                    count = count+1;
                col = [0 0 1; 0 .5 1; 0 1 0;1 .6 0;1 0 0; .5 0 0  ];
                for k = 1:length(types)
                    trials_ktype=(find(trialtypes==types(k)));  
                    temp_data=newrois(trials_ktype,1:length(ts),rois(i));
                    temp_dKappa = dKappa(trials_ktype,:);
                    [all_data,detected] = detect_Ca_events(temp_data,frametime,80);
                    detected_data= all_data(find(detected),:);
                    detected_dKappa = temp_dKappa(find(detected),:);
                    max_dFF=zeros(size(detected_data,1),1);
                    total_dKappa = zeros(size(detected_data,1),1);
                    
                     max_dFF=max(detected_data,[],2);
                    total_dKappa =  sum(detected_dKappa,2);
   

%                      plot(total_dKappa,max_dFF,'Marker','o','color',col(types(k),:),'Markersize',6);
                      scatter(total_dKappa,max_dFF,20,col(types(k),:),'fill'); hold on;
%                     vline([.5 1 1.5 2 2.5],'k-');

                end
%                 axis([0 2000 0 600]);set(gca,'YMinorTick','on','YTick', 0:200:600);
                set(gca,'Xscale','log');

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

function [event_detected_data,detected] = detect_Ca_events(src_data,sampling_time,threshold)
    event_detected_data = zeros(size(src_data,1),floor(3/sampling_time) );
    detected = zeros(size(src_data,1),1);
    for i = 1:size(src_data,1)
       current_trace = src_data(i,:);
       
       if (max(current_trace)> threshold)
           
           [ind,hval] = hist(current_trace);
           F0 = mean(hval(1:4));
           Fmax = max(hval);
           F1 = F0 * ones(round(.5/sampling_time),1);
           F2 = [F0 * ones(round(length(F1)/2)-1,1);Fmax * ones(round(length(F1)/2),1)];
            for j = 1: length(current_trace)
               tp_trace = current_trace(j:j + length(F1)-1);
               lse1 = sum((tp_trace-F1').^2);
               lse2 = sum((tp_trace-F2').^2);
               if(lse2<lse1)
                   detected(i) =1;
                   event_index = j+round(length(F1)/2); 
                   temp = current_trace(max(1,event_index - floor(.5/sampling_time)): min(length(current_trace),event_index + round(2.5/sampling_time)));
                   leading_blank = event_index - floor(.5/sampling_time);
                   offset = (leading_blank <0);
                   event_detected_data (i, (offset*leading_blank*-1)+1 :length(temp)+(offset*leading_blank*-1)) = temp;
                   
                   break
               end
               
               if(j+1>length(current_trace)-length(F1)+1)
                   break
               end
     
            end

       end

    end
    
end

