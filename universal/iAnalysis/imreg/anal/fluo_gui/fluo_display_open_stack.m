%
% Loads a stack into fluo_display - if passed, opens what you pass; otherwise, dialog
%
%  filepath: directory where file resides
%  filename: name of file
%  img_stack: if specified, no loading takes place; this is used ; [] means load 
%  denovo: if 1, then this wipes all processor step images, leaving steps intact but unexecuted.
%          this is used if the user invoes load.
%
function fluo_display_open_stack (filepath, filename, img_stack, denovo)
  global glovars;
	loaded_from_file = 0;

	% denovo?
	if (nargin < 4) ; denovo = 0; end

  % 0) select via dialog ? 
	if (exist('filepath','var') == 0 & exist('filename','var') == 0)
	  disp('No filename specified; invoking UI');
		[filename, filepath]=uigetfile({'*.tif;*.tiff', 'TIFF file (*.tif, *.tiff)'; ...
                     '*.*','All Files (*.*)'},'Select TIFF file', glovars.data_lastpath);
		denovo = 1;
	elseif (exist('filename','var') == 0)
	  disp('Select a file in the specified directory...');
		[filename, filepath]=uigetfile({'*.tif;*.tiff', 'TIFF file (*.tif, *.tiff)'; ...
                     '*.*','All Files (*.*)'},'Select TIFF file', filepath);
		denovo = 1;
	end

	glovars.fluo_display.im_filepath = filepath;
	glovars.fluo_display.im_filename = filename;

	% 1) pull necessary stuff from fluo_control_main
	glovars.fluo_display.nchan = str2num(get(glovars.fluo_control_main.nchan_edit,'String'));
	glovars.fluo_display.usedchan = str2num(get(glovars.fluo_control_main.usedchan_edit,'String'));

  % 2) if valid, load it                 
  fullname = [filepath filesep filename];
	glovars.fluo_display.fluo_display_filepath = filepath;
	if (exist('img_stack','var') == 0)
	  img_stack = [];
	end
	if (length(img_stack) == 0) 
    disp(['Attempting to read ' fullname ' ...']);      

    n_opt.numchannels = glovars.fluo_display.nchan;
		n_opt.channel = glovars.fluo_display.usedchan;

		[fluo_display_im improps] = load_image(fullname, -1, n_opt);
    fluo_display_im_hdr = improps.header;

    % sanity check
		if (length(fluo_display_im) == 0) ; return ; end

		loaded_from_file = 1;
	else
	  fluo_display_im = img_stack;
		fluo_display_im_hdr = 'Created by fluo_display_open_stack.';
	end

  % assign image parameters
	glovars.fluo_display.display_im = fluo_display_im;
	glovars.fluo_display.display_im_hdr = fluo_display_im_hdr;
	glovars.fluo_display.display_im_frame = 1;
  glovars.fluo_display.colormap_min = max(0,min(min(min(fluo_display_im))));
  glovars.fluo_display.colormap_max = max(max(max(fluo_display_im)));
  glovars.fluo_display.display_im_height = size(fluo_display_im,1);
  glovars.fluo_display.display_im_width = size(fluo_display_im,2);
  glovars.fluo_display.display_im_nframes = size(fluo_display_im,3);

	if (loaded_from_file)
	  glovars.fluo_display.nchan = improps.nchans;
	end

  if (glovars.fluo_display.colormap_max > glovars.fluo_display.max_pixel_value)
	  disp(['fluo_display_open_stack.m::max color value exceeds ' num2str(glovars.fluo_display.max_pixel_value) ' ; truncating.']);
		glovars.fluo_display.colormap_max = glovars.fluo_display.max_pixel_value;
	end
  % construct the x and y matrices
	w = size(fluo_display_im,1);
	h = size(fluo_display_im,2);
	glovars.fluo_display.display_im_x_matrix = zeros(w,h);
	glovars.fluo_display.display_im_y_matrix = zeros(w,h);
  for x=1:w
	  glovars.fluo_display.display_im_y_matrix(x,:) = x*ones(1,h);
	end
  for y=1:h
	  glovars.fluo_display.display_im_x_matrix(:,y) = y*ones(w,1);
	end

  % 3) slider

  % 4) cleanup / display
  fluo_display_update_display();

	% 5) If first load, add to processor list
	if (length(glovars.processor_step)== 0 | denovo)
		glovars.processor_step(1).name = 'Base Image';
		glovars.processor_step(1).uid = 0; % always for base image!
		glovars.processor_step(1).im_path = sprintf('%s',glovars.fluo_display.im_filepath);
		glovars.processor_step(1).im_fname = sprintf('%s',glovars.fluo_display.im_filename);
		glovars.processor_step(1).gui_handle = [];
		glovars.processor_step(1).processor.type_id = 0;
		glovars.processor_step(1).processor.name = 'Base Image';
		glovars.processor_step(1).processor.gui_present = 0;
		glovars.processor_step(1).processor.image_output = 1;
		glovars.current_processor_step = 1;

		% set data last path
    glovars.data_lastpath =  sprintf('%s',glovars.fluo_display.im_filepath);
 
    % clear any other steps ...
		for s=2:length(glovars.processor_step)
			glovars.processor_step(s).im_path = '';
			glovars.processor_step(s).im_fname = '';
		end

		fluo_control_main_update_processing_steps_list();
  end

	% 6) roi update
  fluo_roi_control_update_rois();
