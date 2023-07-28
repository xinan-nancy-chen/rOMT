function [path_par,data_par,model_par] = getParams(data_name)
% Created by Rena Elkin on 12/20/2016.
% Last modified: 12/20/2016.
%                01/17/2017.
%                03/20/2021 by Xinan Chen
% This function returns parameters corresponding to the input 'data_name'
%   Input:  'data_name'  := data id specification (string)
%   Output: 'path_par' := structure containing path parameters
%           'data_par' := structure containing data (img) parameters
%           'model_par' := structure containing model parameters

if nargin<1
    data_name='C294';
end
path_par.out_dir_pre = '../rOMT';

switch data_name
    %DOWNSAMPLED VERSION
    case 'C294'%%WT
        path_par.data_tag = 'C294';
        path_par.data_dir = './data/12_MONTH_DATA/psnrv/WT/C294/C294_031318A/';
        path_par.data_name = 'psnrv_C294_031318A_DOTA37_30ul_E';
        path_par.data_name_post = '';
        path_par.extension = '.nii';
        
        path_par.max_dpsnrv = './data/12_MONTH_DATA/MAXpsnrv/C294_031318A_psnrv_max.nii';
        path_par.anato = './data/12_MONTH_DATA/psnrv/WT/C294/pbase_snrv_C294_031318A_DOTA37_30ul_E53.nii';
        
        path_par.sp_mask_opts(1).name = 'brain';
        path_par.sp_mask_opts(1).path = './data/12_MONTH_DATA/12months_mask_brainCSF/C294.nii';
        
        data_par.dataset_name = 'CAA';
        data_par.name = 'C294';
        data_par.mask_number = 1;
        switch data_par.mask_number
            case 1%head mask
                path_par.data_mask_path = './data/12_MONTH_DATA/masksforOMT/WT/C294_MASK.nii';
                data_par.x_range = 20:80;%22:78;
                data_par.y_range = 1:106;%1:106;
                data_par.z_range = 39:85;%41:83;
        end
        
        data_par.domain_size = [100,106,100];
        
        data_par.n1=length(data_par.x_range);
        data_par.n2=length(data_par.y_range);
        data_par.n3=length(data_par.z_range);
        
        data_par.do_time_interp = 0;
        
        switch data_par.do_time_interp
            case 1
                tinterp_str = 'gauss';
            case 0
                tinterp_str = 'none';
        end
        data_par.do_resize = 1;
        data_par.size_factor = 0.5;
        
        data_par.data_index_E = 19:53;
        data_par.inflow_start_time = 22;
        data_par.inflow_stop_time = 37;
        data_par.true_size=round(data_par.size_factor*[length(data_par.x_range),length(data_par.y_range),length(data_par.z_range)]);
        
        data_par.up_thresh=30000;
        data_par.low_thresh=100;
        data_par.redistribute = 0;
        data_par.normalize = 0;
        data_par.smooth = 1;
        data_par.reinitR = 0; %0 if do consecutively and 1 if reinitialize rho
        data_par.dilate = 3;
        
        data_par.first_time = data_par.data_index_E(13);
        data_par.time_jump = 3;
        data_par.last_time = data_par.data_index_E(31);
        
        model_par.sigma = 2e-3;
        model_par.sig_str = '2e3';
        model_par.dt = 0.4;
        model_par.nt = 10;
        model_par.gamma = 0.008;
        model_par.beta  = 0.0001;
        model_par.niter_pcg = 20;
        
        model_par.dTri = 1;%3;%1 := 'closed', 3:= 'open'
        model_par.add_source = 0;
        path_par.version = sprintf('diff_%s_mask_%d_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_uini_0_rini_%s_beta_%5.4f_R_gamma_%4.3f_correctHu_dtri%d_tinterp%d_rmin%d_rnorm%d_rsmooth%d_rreinit%d_source%d_dilate%d_pcg%d',...
            model_par.sig_str,data_par.mask_number,data_par.time_jump,model_par.dt,model_par.nt,data_par.first_time,data_par.last_time,tinterp_str,model_par.beta,model_par.gamma,model_par.dTri,data_par.do_time_interp,data_par.redistribute,data_par.normalize,data_par.smooth,data_par.reinitR,model_par.add_source,data_par.dilate,model_par.niter_pcg);
        path_par.out_dir=sprintf('%s/test_results/%s/%s',path_par.out_dir_pre,path_par.data_tag,path_par.version);
        %%   
    %ORIGINAL VERSION
    %{
    case 'C294'%%WT
        path_par.data_tag = 'C294';
        path_par.data_dir = './data/12_MONTH_DATA/psnrv/WT/C294/C294_031318A/';
        path_par.data_name = 'psnrv_C294_031318A_DOTA37_30ul_E';
        path_par.data_name_post = '';
        path_par.extension = '.nii';
        
        path_par.max_dpsnrv = './data/12_MONTH_DATA/MAXpsnrv/C294_031318A_psnrv_max.nii';
        path_par.anato = './data/12_MONTH_DATA/psnrv/WT/C294/pbase_snrv_C294_031318A_DOTA37_30ul_E53.nii';
        
        path_par.sp_mask_opts(1).name = 'brain';
        path_par.sp_mask_opts(1).path = './data/12_MONTH_DATA/12months_mask_brainCSF/C294.nii';
        
        data_par.dataset_name = 'CAA';
        data_par.name = 'C294';
        data_par.mask_number = 1;
        switch data_par.mask_number
            case 1%head mask
                path_par.data_mask_path = './data/12_MONTH_DATA/masksforOMT/WT/C294_MASK.nii';
                data_par.x_range = 20:80;%22:78;
                data_par.y_range = 1:106;%1:106;
                data_par.z_range = 39:85;%41:83;
        end
        
        data_par.domain_size = [100,106,100];
        
        data_par.n1=length(data_par.x_range);
        data_par.n2=length(data_par.y_range);
        data_par.n3=length(data_par.z_range);
        
        data_par.do_time_interp = 0;
        
        switch data_par.do_time_interp
            case 1
                tinterp_str = 'gauss';
            case 0
                tinterp_str = 'none';
        end
        data_par.do_resize = 0;
        data_par.size_factor = 1;
        
        data_par.data_index_E = 19:53;
        data_par.inflow_start_time = 22;
        data_par.inflow_stop_time = 37;
        data_par.true_size=round(data_par.size_factor*[length(data_par.x_range),length(data_par.y_range),length(data_par.z_range)]);
        
        data_par.up_thresh=30000;
        data_par.low_thresh=100;
        data_par.redistribute = 0;
        data_par.normalize = 0;
        data_par.smooth = 1;
        data_par.reinitR = 0; %0 if do consecutively and 1 if reinitialize rho
        data_par.dilate = 3;
        
        data_par.first_time = data_par.data_index_E(13);
        data_par.time_jump = 2;
        data_par.last_time = data_par.data_index_E(33);
        
        model_par.sigma = 2e-3;
        model_par.sig_str = '2e3';
        model_par.dt = 0.4;
        model_par.nt = 10;
        model_par.gamma = 0.008;
        model_par.beta  = 0.0001;
        model_par.niter_pcg = 60;
        
        model_par.dTri = 1;%3;%1 := 'closed', 3:= 'open'
        model_par.add_source = 0;
        path_par.version = sprintf('diff_%s_mask_%d_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_uini_0_rini_%s_beta_%5.4f_R_gamma_%4.3f_correctHu_dtri%d_tinterp%d_rmin%d_rnorm%d_rsmooth%d_rreinit%d_source%d_dilate%d_pcg%d',...
            model_par.sig_str,data_par.mask_number,data_par.time_jump,model_par.dt,model_par.nt,data_par.first_time,data_par.last_time,tinterp_str,model_par.beta,model_par.gamma,model_par.dTri,data_par.do_time_interp,data_par.redistribute,data_par.normalize,data_par.smooth,data_par.reinitR,model_par.add_source,data_par.dilate,model_par.niter_pcg);
        path_par.out_dir=sprintf('%s/results/CAA12m/%s/%s',path_par.out_dir_pre,path_par.data_tag,path_par.version);
        %%
    %}
    otherwise
        disp('data-name unrecognized')
end
