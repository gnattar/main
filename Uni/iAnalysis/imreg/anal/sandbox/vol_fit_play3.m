
%% iterate over different angles, z's
im_list = {'an160508_2012_02_10_based_2012_02_10_fov_01002','an160508_2012_02_10_based_2012_02_10_fov_01003','an160508_2012_02_10_based_2012_02_10_fov_01004', ...
       'an160508_2012_02_10_based_2012_02_10_fov_02002','an160508_2012_02_10_based_2012_02_10_fov_02003','an160508_2012_02_10_based_2012_02_10_fov_02004', ...
       'an160508_2012_02_10_based_2012_02_10_fov_03002','an160508_2012_02_10_based_2012_02_10_fov_03003','an160508_2012_02_10_based_2012_02_10_fov_03004', ...
       'an160508_2012_02_10_based_2012_02_10_fov_04002','an160508_2012_02_10_based_2012_02_10_fov_04003','an160508_2012_02_10_based_2012_02_10_fov_04004', ...
       'an160508_2012_02_10_based_2012_02_10_fov_05002','an160508_2012_02_10_based_2012_02_10_fov_05003','an160508_2012_02_10_based_2012_02_10_fov_05004', ...
       'an160508_2012_02_10_based_2012_02_10_fov_06002','an160508_2012_02_10_based_2012_02_10_fov_06003','an160508_2012_02_10_based_2012_02_10_fov_06004', ...
       'an160508_2012_02_10_based_2012_02_10_fov_07002','an160508_2012_02_10_based_2012_02_10_fov_07003','an160508_2012_02_10_based_2012_02_10_fov_07004', ...
       'an160508_2012_02_10_based_2012_02_10_fov_08002','an160508_2012_02_10_based_2012_02_10_fov_08003','an160508_2012_02_10_based_2012_02_10_fov_08004', ...
       'an160508_2012_02_10_based_2012_02_10_fov_09002','an160508_2012_02_10_based_2012_02_10_fov_09003','an160508_2012_02_10_based_2012_02_10_fov_09004', ...
       'an160508_2012_02_10_based_2012_02_10_fov_10002','an160508_2012_02_10_based_2012_02_10_fov_10003','an160508_2012_02_10_based_2012_02_10_fov_10004', ...
       'an160508_2012_02_10_based_2012_02_10_fov_11002','an160508_2012_02_10_based_2012_02_10_fov_11003','an160508_2012_02_10_based_2012_02_10_fov_11004'};

for i=1:length(im_list)
  im_list{i} = [pwd filesep 'rois' filesep im_list{i}];
end

cd /media/an160508b;
if (~exist('vim'))
	vim = single(load_image('an160508_2012_02_08_dz1_stack_0to611_pz250_start_15_001'));
end
whos('vim');

guess_z = 80:15:80+(15*(length(im_list)-1));
out_path_root = '/media/an160508b/volume_fit';
fit_all_planes_to_vol_warpfield(vim, im_list, guess_z, out_path_root);

