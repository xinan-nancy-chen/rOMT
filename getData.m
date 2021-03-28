function im_out = getData(data_tag,tind,tinterp)
%this function loads the data images with no processing or interpolation
%last modified: 01/17/2017
%               01/23/2017
%               01/26/2017

addpath('../../Visualization/MATLAB MRI display/NIfTI_analyze/');

if nargin < 1
    data_tag = 'C294';
    tind = 16;
    tinterp = 'none';
elseif nargin < 2
    tind = 16;
    tinterp = 'none';
elseif nargin < 3
    tinterp = 'none';
end

[path_par,data_par]=getParams(data_tag);
if isfield(data_par,'dataset_name')
    switch data_par.dataset_name
        case {'LN'} % 4D data
            D=load_untouch_nii(sprintf('%s%s%s',path_par.data_dir,path_par.data_name,path_par.extension));
            im = squeeze(D.img(:,:,:,tind));
        otherwise
            D=load_untouch_nii(sprintf('%s%s%02d%s%s',path_par.data_dir,path_par.data_name,tind,path_par.data_name_post,path_par.extension));
            im=D.img;
    end
            
else
    switch path_par.data_tag
        case {'tumor'}
            D=load_untouch_nii(sprintf('%s%s',path_par.data_dir,path_par.extension));
            im = squeeze(D.img(:,:,:,tind));
        otherwise
            D=load_untouch_nii(sprintf('%s%s%02d%s%s',path_par.data_dir,path_par.data_name,tind,path_par.data_name_post,path_par.extension));
            im=D.img;
    end
end

if data_par.do_time_interp
    time_range = data_par.first_time:data_par.last_time;
    switch tinterp
        case 'none'
            a=0;
        case {'Gaussian', 'Gauss', 'gaussian', 'gauss'}
            weights = [ 0.6065    0.1353    0.0111 ]; %from sigma = 1 gaussian weights: [0.0111    0.1353    0.6065    1.0000    0.6065    0.1353    0.0111]
            tot_weight = 1;
            
            for ind_weight = 1:min([length(weights), tind-1-(time_range(1)), (time_range(end))-tind])
                D1_add = load_untouch_nii(sprintf('%s%s%d%s%s',path_par.data_dir,path_par.data_name,tind+ind_weight,path_par.data_name_post,path_par.extension), '');
                im = im + weights(ind_weight)*D1_add.img;
                
                D1_add = load_untouch_nii(sprintf('%s%s%d%s%s',path_par.data_dir,path_par.data_name,tind-ind_weight,path_par.data_name_post,path_par.extension), '');
                im = im + weights(ind_weight)*D1_add.img;
                
                tot_weight = tot_weight + 2*weights(ind_weight);
            end
            im = im/tot_weight;
    end
end

im=double(im(data_par.x_range,data_par.y_range,data_par.z_range));
im(im<0)=0;
im(im>data_par.up_thresh) = data_par.up_thresh;

if data_par.do_resize
    S_in = size(im);
    S_out = round(S_in*data_par.size_factor);
    im = resizeMatrix(im,S_out,'linear');%resizeMatrix(im,S_out,'cubic');
end


if data_par.redistribute 
    min_val=data_par.up_thresh*0.1;
    im = im+min_val;
    im = im./(data_par.up_thresh + min_val);
end

if data_par.normalize
    m = nii2mat(path_par.data_mask_path,data_par.x_range,data_par.y_range,data_par.z_range);
    imM = im; imM(~m) = 0;
    im = im.*data_par.up_thresh./sum(imM(:));
end


im_out=im;

end
