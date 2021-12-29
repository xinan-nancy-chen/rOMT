function gladpar = getGLADpar(tag,type)
% Created by Xinan Chen on 03/24/2021.
% This function returns parameters corresponding to the input 'data_name' and analysis type
%   Input:  type  := 's' if run for speed and Peclet maps; 'v' if run for pathlines w./w.o. endowed values
%   Output: 'glad_par' := structure containing parameters for GLAD analysis

%%
[ppar,dpar,mpar] = getParams(tag);
gladpar.do_spMsk = 0; %1 if use alternative mask for selecting pathline start points
gladpar.spMsk_name = 'olf_nomain'; %name of mask to use for selecting start points

% find index of mask with spMsk_name
% get path of spMsk
% check if spMsk is the same as the msk, if so -> set do_spMsk = 0
if gladpar.do_spMsk && isfield(ppar,'sp_mask_opts')
    sp_msk_path_opts = {ppar.sp_mask_opts(:).name};
    gladpar.spMsk_ind = find(strcmp(gladpar.spMsk_name,sp_msk_path_opts));
    if isempty(gladpar.spMsk_ind)
        warning(fprintf('getLPPA.m: unable to find path coresponding to mask named %s\n>>> Default mask being used and ''do_spMsk'' set to 0\n',gladpar.spMsk_name));
        gladpar.do_spMsk = 0; 
        clear gladpar.spMsk_name gladpar.spMsk_ind
    end        
else
    gladpar.do_spMsk = 0; 
    clear gladpar.spMsk_name gladpar.spMsk_ind
end
if gladpar.do_spMsk
    gladpar.spMSK_str = sprintf('altSPmsk_%s_opt%d',gladpar.spMsk_name,gladpar.spMsk_ind);
else
    gladpar.spMSK_str = 'altSPmsk0';
end
    
gladpar.do_spErodeMsk = 0; %0 if don't erode SPmask to find start points, N if erode by N points
gladpar.do_spDilateMsk = 0; %0 if don't dilate SPmask for start points

if gladpar.do_spErodeMsk>0 && gladpar.do_spDilateMsk>0
    warning('getGLADpar.m: spMask cannot be both eroded and dilated. spMask will be eroded (no dilation) by default');
    gladpar.do_spDilateMsk = 0;
end

gladpar.do_masked = 1;%1 if ensure everything is inside mask before saving to .nii

%% set parameters, start points and directory
switch dpar.dataset_name
    case 'CAA'
        switch type
            case 's'
                gladpar.do_sp = 1; %1 if use sp from max(dpsnrv)>sp_thresh, 0 o/w
                gladpar.sp_thresh = 12;
                
                gladpar.mdt = 10;
                gladpar.sl_tol = 2; %threshold for minimum Euclidean length between initial and final streamline points                
                gladpar.spType = 'ordered'; %'uniform' (default: 'ordered'); %'ordered': every fs-th possible sp is used (default)sp distribution; 'uniform': sp chosen according to uniform distribution
                switch gladpar.spType
                    case 'ordered'
                        gladpar.fs = 5;
                        gladpar.spSTR = sprintf('spDISTordered_fs%d',gladpar.fs);
                        gladpar.spSTRtitle = sprintf('spDISTordered-fs%d',gladpar.fs);
                    case 'uniform'
                        gladpar.spPerc = 40;
                        gladpar.spSTR = sprintf('spDISTunif_spPerc%d',gladpar.spPerc);
                        gladpar.spSTRtitle = sprintf('spDISTunif-spPerc%d',gladpar.spPerc);
                end
                gladpar.nEulStep = 1; %number of Eulerian steps to be taken

                % QuickBundle parameters 
                gladpar.do_centroid = 0;
                gladpar.qbthresh = 4;%5;%QuickBundles threshold distance
                gladpar.clus_tol = 12; %threshold for min # of sl's needed in a cluster
                gladpar.clus_cutoff = 600; %pick up to # of largest clusters
                gladpar.nbp = 124; %# of points that QuickBundles will resample streamline to
                gladpar.metric_str = 'AveragePointwiseEuclideanMetric';
                switch gladpar.metric_str
                    case 'AveragePointwiseEuclideanMetric'
                        gladpar.ms = 'AvgPwEuc';
                        gladpar.feature_str = sprintf('ResampleFeature(nb_points=%d)',gladpar.nbp);
                    case 'CosineMetric'
                        gladpar.ms = 'Cos';
                        gladpar.feature_str = 'VectorOfEndpointsFeature()';
                end
                
            case 'v'
                gladpar.do_sp = 0; %1 if use sp from max(dpsnrv)>sp_thresh, 0 o/w
                gladpar.sp_thresh = 12;
                
                gladpar.mdt = 10;
                gladpar.sl_tol = 2; %threshold for minimum Euclidean length between initial and final streamline points                
                gladpar.spType = 'ordered'; %'uniform' (default: 'ordered'); %'ordered': every fs-th possible sp is used (default)sp distribution; 'uniform': sp chosen according to uniform distribution
                switch gladpar.spType
                    case 'ordered'
                        gladpar.fs = 1;
                        gladpar.spSTR = sprintf('spDISTordered_fs%d',gladpar.fs);
                        gladpar.spSTRtitle = sprintf('spDISTordered-fs%d',gladpar.fs);
                    case 'uniform'
                        gladpar.spPerc = 40;
                        gladpar.spSTR = sprintf('spDISTunif_spPerc%d',gladpar.spPerc);
                        gladpar.spSTRtitle = sprintf('spDISTunif-spPerc%d',gladpar.spPerc);
                end
                gladpar.nEulStep = 1; %number of Eulerian steps to be taken
                
                % NCA parameters
                gladpar.radius = 2; % redius to dialate
                gladpar.NNnum_tol = 20; 
                gladpar.stdcut = 0.5; 
                gladpar.meancut = [0.5,1]; 
                gladpar.Npcut = 1; 
                gladpar.Avepathlcut = 8; 
                %gladpar.AveMaxpathscut = 0.02;
                gladpar.threshstr = sprintf('NNnum_tol_%d_||_stdcut_%.2f_wmeancut_[%.2f,%.2f]_Npcut_%d_Avepathlcut_%.1f',gladpar.NNnum_tol,gladpar.stdcut,gladpar.meancut(1),gladpar.meancut(2),gladpar.Npcut,gladpar.Avepathlcut);
                
                gladpar.maskdilate = 3;
                gladpar.radius2 = 3; % redius to dialate
                gladpar.NNnum_tol2 = 10; 
                gladpar.stdcut2 = 0.4; 
                gladpar.meancut2 = [0.6,1]; 
                gladpar.Npcut2 = 1; 
                gladpar.Avepathlcut2 = 10; 
                gladpar.AveMaxpathscut2 = 0.02; 
                gladpar.pathlcut2 = 10;%20;
                gladpar.threshstr2 = sprintf('NNnum_tol_%d_||_stdcut_%.2f_wmeancut_[%.2f,%.2f]_Npcut_%d_Avepathlcut_%.1f|pathlcut_%.1f',gladpar.NNnum_tol2,gladpar.stdcut2,gladpar.meancut2(1),gladpar.meancut2(2),gladpar.Npcut2,gladpar.Avepathlcut2,gladpar.pathlcut2);
                
                
                % vis parameters
                gladpar.jp = 5;
                
                
                
            otherwise
                fprintf('getGLADpar: non-applicable type!');
                return
        end
    %{  
    case 'LN'
        switch type
            case 's'
                gladpar.do_sp = 0; %1 if use sp from max(dpsnrv)>sp_thresh, 0 o/w
                gladpar.sp_thresh = 12;
                
                gladpar.mdt = 10;
                gladpar.sl_tol = 2; %threshold for minimum Euclidean length between initial and final streamline points                
                gladpar.spType = 'ordered'; %'uniform' (default: 'ordered'); %'ordered': every fs-th possible sp is used (default)sp distribution; 'uniform': sp chosen according to uniform distribution
                switch gladpar.spType
                    case 'ordered'
                        gladpar.fs = 5;
                        gladpar.spSTR = sprintf('spDISTordered_fs%d',gladpar.fs);
                        gladpar.spSTRtitle = sprintf('spDISTordered-fs%d',gladpar.fs);
                    case 'uniform'
                        gladpar.spPerc = 40;
                        gladpar.spSTR = sprintf('spDISTunif_spPerc%d',gladpar.spPerc);
                        gladpar.spSTRtitle = sprintf('spDISTunif-spPerc%d',gladpar.spPerc);
                end
                gladpar.nEulStep = 1; %number of Eulerian steps to be taken

                % QuickBundle parameters 
                gladpar.do_centroid = 0;
                gladpar.qbthresh = 4;%5;%QuickBundles threshold distance
                gladpar.clus_tol = 12; %threshold for min # of sl's needed in a cluster
                gladpar.clus_cutoff = 600; %pick up to # of largest clusters
                gladpar.nbp = 124; %# of points that QuickBundles will resample streamline to
                gladpar.metric_str = 'AveragePointwiseEuclideanMetric';
                switch gladpar.metric_str
                    case 'AveragePointwiseEuclideanMetric'
                        gladpar.ms = 'AvgPwEuc';
                        gladpar.feature_str = sprintf('ResampleFeature(nb_points=%d)',gladpar.nbp);
                    case 'CosineMetric'
                        gladpar.ms = 'Cos';
                        gladpar.feature_str = 'VectorOfEndpointsFeature()';
                end
                
            case 'v'
                gladpar.do_sp = 0; %1 if use sp from max(dpsnrv)>sp_thresh, 0 o/w
                gladpar.sp_thresh = 12;
                
                gladpar.mdt = 10;
                gladpar.sl_tol = 1; %threshold for minimum Euclidean length between initial and final streamline points                
                gladpar.spType = 'ordered'; %'uniform' (default: 'ordered'); %'ordered': every fs-th possible sp is used (default)sp distribution; 'uniform': sp chosen according to uniform distribution
                switch gladpar.spType
                    case 'ordered'
                        gladpar.fs = 1;
                        gladpar.spSTR = sprintf('spDISTordered_fs%d',gladpar.fs);
                        gladpar.spSTRtitle = sprintf('spDISTordered-fs%d',gladpar.fs);
                    case 'uniform'
                        gladpar.spPerc = 40;
                        gladpar.spSTR = sprintf('spDISTunif_spPerc%d',gladpar.spPerc);
                        gladpar.spSTRtitle = sprintf('spDISTunif-spPerc%d',gladpar.spPerc);
                end
                gladpar.nEulStep = 1; %number of Eulerian steps to be taken
                
                % NCA parameters
                gladpar.radius = 2; % redius to dialate
                gladpar.NNnum_tol = 20; 
                gladpar.stdcut = 0.5; 
                gladpar.meancut = [0.5,1]; 
                gladpar.Npcut = 1; 
                gladpar.Avepathlcut = 8; 
                %gladpar.AveMaxpathscut = 0.02;
                gladpar.threshstr = sprintf('NNnum_tol_%d_||_stdcut_%.2f_wmeancut_[%.2f,%.2f]_Npcut_%d_Avepathlcut_%.1f',gladpar.NNnum_tol,gladpar.stdcut,gladpar.meancut(1),gladpar.meancut(2),gladpar.Npcut,gladpar.Avepathlcut);
                
                gladpar.maskdilate = 3;
                gladpar.radius2 = 3; % redius to dialate
                gladpar.NNnum_tol2 = 10; 
                gladpar.stdcut2 = 0.4; 
                gladpar.meancut2 = [0.6,1]; 
                gladpar.Npcut2 = 1; 
                gladpar.Avepathlcut2 = 10; 
                gladpar.AveMaxpathscut2 = 0.02; 
                gladpar.pathlcut2 = 20;
                gladpar.threshstr2 = sprintf('NNnum_tol_%d_||_stdcut_%.2f_wmeancut_[%.2f,%.2f]_Npcut_%d_Avepathlcut_%.1f|pathlcut_%.1f',gladpar.NNnum_tol2,gladpar.stdcut2,gladpar.meancut2(1),gladpar.meancut2(2),gladpar.Npcut2,gladpar.Avepathlcut2,gladpar.pathlcut2);
                
                
                % vis parameters
                gladpar.jp = 5;
                
            otherwise
                fprintf('getGLADpar: non-applicable type!');
                return
        end
    %}    
    otherwise
        fprintf('getGLADpar: non-applicable dataset_name!');
        return

end


gladpar.flw_type = 'vel';%'flw';
gladpar.pln = 2; %minimum number of unique points needed to be considered a pathline
gladpar.RD = 'R'; %'D' if use data; 'R' if use interpolated img
gladpar.XT = 'T'; %'X' if interp fixed spatial dist, 'T' if use time step;
gladpar.minIm0 = 0; %1 if set img = img - min(img), 0 otherwise

gladpar.cutoff_str = 'min';
switch gladpar.cutoff_str
    case 'min'
        gladpar.thresholds.conc = 0.0001;
        gladpar.thresholds.front = 0;
        gladpar.thresholds.flow = 0;
        gladpar.thresholds.speed = 0.0001;
    case 'max'
        gladpar.thresholds.conc = .0001;
        gladpar.thresholds.front = 0;
        gladpar.thresholds.flow = 0;
        gladpar.thresholds.speed = 0;
    case 'mean'
        gladpar.thresholds.conc = .00001;
        gladpar.thresholds.front = 0;
        gladpar.thresholds.flow = 0.0001;
        gladpar.thresholds.speed = 0.0001;
end



end