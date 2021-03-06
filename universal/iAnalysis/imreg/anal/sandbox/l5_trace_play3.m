%
%
% start_z = 368 ; end_z = 619 ;
% tim = mean(vim(:,:,start_z-5:start_z+5),3);
% figure ; imshow(tim,[0 1000]);
% rA = roi.roiArray();
% rA.masterImage = tim;
% rA.startGui;
%
% select it . . . 
%
% l5_trace_play2(vim, rA, start_z, end_z)
%
function l5_trace_play2(rO, rId, vim)
% now try to follow it thru the stack . . .

  % Prep work
	rA = rO;
	rIdx = find(rO.roiIds == rId);
	nPixels = length(rO.rois{rIdx}.indices);
	allRs = {};

  % map to volume
  [col row] = get_xy_from_idx(rO.rois{rIdx}.indices, size(rO.refStackMap.X));
	X = 0*length(col);
	Y = 0*length(col);
	Z = 0*length(col);
	for i=1:length(col)
  	X(i) = rO.refStackMap.X(row(i),col(i));
  	Y(i) = rO.refStackMap.Y(row(i),col(i));
  	Z(i) = rO.refStackMap.Z(row(i),col(i));
	end
	X = round(X);
	Y = round(Y);
	Z = round(Z);

  % load volume
%	vim = load(assign_root_data_path(rO.refStackRelPath));
  end_z = size(vim,3);

  % base image yo
	start_z = round(mean(Z));
	tRoi = roi.roi(1, [], get_idx_from_xy(X,Y, size(rO.refStackMap.X)), [1 1 0], size(rO.refStackMap), []);
	tRoi.imageBounds = size(rO.refStackMap.X);
	tRoi.indices = get_idx_from_xy(X,Y, size(rO.refStackMap.X));
	tRoi.corners =tRoi.computeBoundingPoly(); 
	tRoi.assignBorderIndices();

	rA = roi.roiArray();
	rA.masterImage = vim(:,:,start_z);
	rA.workingImage = vim(:,:,start_z);
	rA.addRoi(tRoi);
	%rA.startGui();
%	return

% gradient autodet settings
gparams.radRange =  [3 10];
gparams.edgeSign = -1;
gparams.postDilateBy = 0;
gparams.dPosMax = 2;
gparams.postFracRemove = 0;

% START AT THE BOTTOM WORK UP
  end_z = 300;
	fin_z = end_z;
	close all;

	ax = axes;
	nPixels = [];
	for i=start_z:-1:end_z
    try 
			nim = mean(vim(:,:,i-5:i+5),3);

			rB = roi.roiArray();
			rB.masterImage = nim;
			rB.workingImage = nim;

			% 1) find best image match in next step based on LAST step
      fm_params(1).value = [];
			fm_params(2).value = [];
			fm_params(3).value = [];
      fm_params(4).value = 5/512;
			fm_params(5).value = [];
			corrvals(i) = roi.roiArray.findMatchingRoisInNewImage(rA, rB, fm_params);

centerPoint0(1) = mean(rA.rois{1}.corners(1,:));
centerPoint0(2) = mean(rA.rois{1}.corners(2,:));
disp(' ');
disp('=================================================================');
disp(['x0: ' num2str(centerPoint0(1)) ' y0: ' num2str(centerPoint0(2))]);

      % 2) 
			nPixels = max(10,length(rA.rois{1}.indices)-5);
			valVec = reshape(rA.workingImage(rA.rois{1}.indices),[],1);
			cutoff = quantile(valVec, 0.2);

			% No growth!
			rB.dilateRoiBorders(1,1);
			rB.fillToBorder(1);
			[irr sidx] = sort(rB.workingImage(rB.rois{1}.indices), 'descend');
			rB.rois{1}.indices = rB.rois{1}.indices(sidx(1:min(nPixels,length(rB.rois{1}.indices))));
			rB.removeAllButLargestConnectedComponent(1);
			rB.fitRoiBordersToIndices(1);

centerPoint0(1) = mean(rB.rois{1}.corners(1,:));
centerPoint0(2) = mean(rB.rois{1}.corners(2,:));
disp(['x: ' num2str(centerPoint0(1)) ' y: ' num2str(centerPoint0(2))]);
			
			% wrapup
			rA = rB;
			tim = nim;
			allRs{i} = rB;

			% plot
			if (1)
				oWi = allRs{i}.workingImage;
				nWi = oWi;
				nWi(find(nWi > 700)) = 700;
				allRs{i}.workingImage = nWi;
				allRs{i}.plotImage(0,1,0,[],ax);
				text(30,30, num2str(i), 'Color', [1 1 0]);
				drawnow;
				allRs{i}.workingImage = oWi;
				disp(['done with ' num2str(i)]);
			end
		catch ME
		  disp(ME.getReport);
		  break; 
		end
	end

fin_z = i+1;
length(allRs)
% build the soma so it is nice . . .

% plot
	close all;
%	figure ; plot(corrvals);
aviobj = avifile(['~/Desktop/' num2str(rId) 'trace.avi']);
fig=figure;
	ax = axes;
	for i=start_z:-1:fin_z
	  oWi = allRs{i}.workingImage;
		nWi = oWi;
	  nWi(find(nWi > 700)) = 700;
		allRs{i}.workingImage = nWi;
	  allRs{i}.plotImage(0,1,0,[],ax);
		text(30,30, num2str(i), 'Color', [1 1 0]);
		drawnow;
		F = getframe(fig);
		aviobj = addframe(aviobj,F);
		allRs{i}.workingImage = oWi;
	end
aviobj = close(aviobj);


% TO DO: spin the cell!
% TO DO: trace UPWARDS
% TO DO: multiple @ once?

% eventual workflow:
% 1) select apical ROIs in an imaged plane
% 2) correlate that plane to the volume stack [volume PATH as roiArray field ; roiArray method registerToStack which does either workingImage or masterImage]
% 3) trace the apical ROIs to their somas and record the soma position [roiArray method traceApicalDendriteThruStack; roiArray field to store full cell?]
% 4) for other planes, determine where the apical dendrite is and label accordingly.
