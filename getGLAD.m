%% This is the the code for post-processing rOMT results
% contains GLAD analysis + NCA

addpath('./Inverse','./Sensitivities','./analyzeFlows',genpath('./utilities'))


%% change input here
% make sure in the path directory in getParams.m is what the post-processing you
% would like to run
tag = 'C294';



%%
fprintf('\n =============== Post-processing Starts ===============\n')
fprintf('_________________________________________________________\n\n')
prompt = '\n=== What result do you want to extract? === \nType "s" for speed and Peclet maps, "v" for pathlines w./w.o. endowed values, "b" for both. \nThen press Enter to continue\n';
str = input(prompt,'s');

switch str %'s' if run for speed and Peclet maps; 'v' if run for pathlines w./w.o. endowed values
    case 's'
        types = {'s'};
    case 'v'
        types = {'v'};
    case 'b'
        types = {'s','v'};
    otherwise
        disp('Input unrecognized, please re-run')
        return
end
%%
for l = 1:length(types)
type = types{l};
switch type 
    case 's'
        paper_fig_str = sprintf('set0%02d',1);
    case 'v'
        paper_fig_str = sprintf('set0%02d',2);
    otherwise
        fprintf('getGLAD: non-applicable type!');
        return
end
formatOut = 'mmddyy';
date_str = datestr(now,formatOut);

[ppar,dpar,mpar] = getParams(tag);

% get GLAD parameters
gladpar = getGLADpar(tag,type);

% set usefule variables:
nt = mpar.nt;
n = dpar.true_size;
ti = dpar.first_time;
tf = dpar.last_time;
tj = dpar.time_jump;

%% load mask:
model_mask = load_untouch_nii(ppar.data_mask_path);
model_mask.hdr.dime.datatype = 16;
model_mask.hdr.dime.bitpix = 32;
if dpar.do_resize
    S = round(dpar.domain_size.*dpar.size_factor);
    model_mask.img = zeros(S);
    model_mask.hdr.dime.dim(2:4) = S;
else
    S = dpar.domain_size;
    model_mask.img = zeros(S);
end

mask2 = nii2mat(ppar.data_mask_path,dpar.x_range,dpar.y_range,dpar.z_range);
msk = zeros(size(mask2));
msk(mask2>0) = 1;
clear mask2

if dpar.do_resize
   msk = resizeMatrix(msk,round(dpar.size_factor.*size(msk)),'linear');
   msk(msk~=1) = 0;
end

if gladpar.do_spMsk
    mskSP = nii2mat(ppar.sp_mask_opts(gladpar.spMsk_ind).path,dpar.x_range,dpar.y_range,dpar.z_range);
    if dpar.do_resize
       mskSP = resizeMatrix(mskSP,round(dpar.size_factor.*size(mskSP)),'linear');
       mskSP(mskSP~=1) = 0;
    end
    mskSP(~msk) = 0; %only include part of spMSK that is inside the whole msk
else
    mskSP = msk;
end


if gladpar.do_spErodeMsk > 0
    [x_er,y_er,z_er] = meshgrid(-gladpar.do_spErodeMsk:gladpar.do_spErodeMsk,-gladpar.do_spErodeMsk:gladpar.do_spErodeMsk,-1:1);
    strel_e = (x_er/gladpar.do_spErodeMsk).^2 + (y_er/gladpar.do_spErodeMsk).^2 + (z_er).^2 <= 1;    
    mskSP = imerode(mskSP, strel_e);
elseif gladpar.do_spDilateMsk > 0    
    %first find if original mask was dilated for rOMT run
    mskDIL = msk;
    if dpar.dilate>0
        [xr,yr,zr] = meshgrid(-dpar.dilate:dpar.dilate,-dpar.dilate:dpar.dilate,-dpar.dilate:dpar.dilate);
        strel = (xr/dpar.dilate).^2 + (yr/dpar.dilate).^2 + (zr/dpar.dilate).^2 <= 1;
        mskDIL = imdilate(msk,strel);
    end
    %now dilate spMsk and remove any points outside of mask used for rOMT
    [x_dil,y_dil,z_dil] = meshgrid(-gladpar.do_spDilateMsk:gladpar.do_spDilateMsk,-gladpar.do_spDilateMsk:gladpar.do_spDilateMsk,-gladpar.do_spDilateMsk:gladpar.do_spDilateMsk);
    strel_d = (x_dil/gladpar.do_spDilateMsk).^2 + (y_dil/gladpar.do_spDilateMsk).^2 + (z_dil/gladpar.do_spDilateMsk).^2 <= 1;
    mskSP = imdilate(mskSP, strel_d);
    mskSP(~mskDIL) = 0;    
end
%mskSP(~msk) = 0;  

if gladpar.do_sp
    dpsnrv_max = nii2mat(ppar.max_dpsnrv,dpar.x_range,dpar.y_range,dpar.z_range);
    if dpar.do_resize
       dpsnrv_max = resizeMatrix(dpsnrv_max,round(dpar.size_factor.*size(dpsnrv_max)),'linear');
    end
    mind = find((mskSP>0) & (dpsnrv_max>gladpar.sp_thresh));
    gladpar.do_sp_str = sprintf('_dpsnrv_min_%d',gladpar.sp_thresh);
else
    mind = find(mskSP>0);
    gladpar.do_sp_str = '';
end

sig_str = mpar.sig_str;
outdir = sprintf('LPPA_%s_%s',paper_fig_str,date_str);
switch type
    case 's'
        outdir_long = sprintf('LPPA_%s_mask%d_type_%s_%s_img%s_mdt%d%s_%s_erode%d_dilate%d_nEul%d_cutoffStr_%s_concThresh%5.4f_spdThresh%5.4f_minIm0_%d_%s_slTol%d_clusters_centroid%d_nbp%d_metric%s_qbthresh%d_clusTol%d_clusCutoff%d_diff%s_tj%d_nt%d%s_%s_%s',...
            tag,dpar.mask_number,'speedmap',gladpar.flw_type,gladpar.RD,gladpar.mdt,gladpar.XT,gladpar.spMSK_str,gladpar.do_spErodeMsk,gladpar.do_spDilateMsk,gladpar.nEulStep,gladpar.cutoff_str,gladpar.thresholds.conc,gladpar.thresholds.speed,gladpar.minIm0,gladpar.spSTR,gladpar.sl_tol,gladpar.do_centroid,gladpar.nbp,gladpar.ms,gladpar.qbthresh,gladpar.clus_tol,gladpar.clus_cutoff,sig_str,tj,nt,gladpar.do_sp_str,paper_fig_str,date_str);
    case 'v'
        outdir_long = sprintf('LPPA_%s_mask%d_type_%s_%s_img%s_mdt%d%s_%s_erode%d_dilate%d_nEul%d_cutoffStr_%s_concThresh%5.4f_spdThresh%5.4f_minIm0_%d_%s_slTol%d_%s_%s_diff%s_tj%d_nt%d%s_%s_%s',...
            tag,dpar.mask_number,'vectors',gladpar.flw_type,gladpar.RD,gladpar.mdt,gladpar.XT,gladpar.spMSK_str,gladpar.do_spErodeMsk,gladpar.do_spDilateMsk,gladpar.nEulStep,gladpar.cutoff_str,gladpar.thresholds.conc,gladpar.thresholds.speed,gladpar.minIm0,gladpar.spSTR,gladpar.sl_tol,gladpar.threshstr,gladpar.threshstr2,sig_str,tj,nt,gladpar.do_sp_str,paper_fig_str,date_str);
    otherwise
        fprintf('getGLAD: non-applicable dataset_name!');
        return
end
if ~exist(sprintf('%s/%s',ppar.out_dir,outdir),'dir')
    mkdir(sprintf('%s/%s',ppar.out_dir,outdir))
end

%%
%
fid = fopen(sprintf('%s/%s/%s_record_%s_%s.txt',ppar.out_dir,outdir,tag,paper_fig_str,date_str),'a+');
fprintf(fid,'%s/%s/%s_record_%s_%s\noutdir-long: %s',ppar.out_dir,outdir,tag,paper_fig_str,date_str,outdir_long);
switch type
    case 's'
        title_str = sprintf('============= initiating...\n\nLagrangian-Pathline (%s data, mask = %d, affSmooth = %d, dilate = %d), \nanalysis type = %s\n \nflw_type = %s, img = %s, mdt = %d(%s), %s, spErode = %d, spDilate = %d, nEulStep= %d, \ncutoffStr = %s, concThresh = %5.4f, spdThresh = %5.4f, minIm0 = %d, %s, slTol = %d, \nclusters, centroid = %d, nbp = %d, metric = %s, qbthresh = %d, clusTol = %d, clusCutoff = %d, \ndiff = %s, tj = %d, nt= %d%s_%s_%s\n\n',...
            tag,dpar.mask_number,dpar.smooth,dpar.dilate,'speedmap',gladpar.flw_type,gladpar.RD,gladpar.mdt,gladpar.XT,gladpar.spMSK_str,gladpar.do_spErodeMsk,gladpar.do_spDilateMsk,gladpar.nEulStep,gladpar.cutoff_str,gladpar.thresholds.conc,gladpar.thresholds.speed,gladpar.minIm0,gladpar.spSTRtitle,gladpar.sl_tol,gladpar.do_centroid,gladpar.nbp,gladpar.ms,gladpar.qbthresh,gladpar.clus_tol,gladpar.clus_cutoff,sig_str,tj,nt,gladpar.do_sp_str,paper_fig_str,date_str);
    case 'v'
        title_str = sprintf('============= initiating...\n\nLagrangian-Pathline (%s data, mask = %d, affSmooth = %d, dilate = %d), \nanalysis type = %s\n \nflw_type = %s, img = %s, mdt = %d(%s), %s, spErode = %d, spDilate = %d, nEulStep= %d, \ncutoffStr = %s, concThresh = %5.4f, spdThresh = %5.4f, minIm0 = %d, %s, slTol = %d, \n%s\n%s\ndiff = %s, tj = %d, nt = %d%s_%s_%s\n\n',...
            tag,dpar.mask_number,dpar.smooth,dpar.dilate,'vectors',gladpar.flw_type,gladpar.RD,gladpar.mdt,gladpar.XT,gladpar.spMSK_str,gladpar.do_spErodeMsk,gladpar.do_spDilateMsk,gladpar.nEulStep,gladpar.cutoff_str,gladpar.thresholds.conc,gladpar.thresholds.speed,gladpar.minIm0,gladpar.spSTRtitle,gladpar.sl_tol,gladpar.threshstr,gladpar.threshstr2,sig_str,tj,nt,gladpar.do_sp_str,paper_fig_str,date_str);
    otherwise
        fprintf('getGLAD: non-applicable dataset_name!');
        return
end
fprintf(title_str)
%%

[x, y, z] = meshgrid(1:n(2), 1:n(1), 1:n(3));
[syind,sxind,szind] = ind2sub(n,mind); %find indices of all voxels inside the ROI-SP
switch gladpar.spType
    case 'ordered'
        sy = syind(1:gladpar.fs:end);
        sx = sxind(1:gladpar.fs:end);
        sz = szind(1:gladpar.fs:end);
    case 'uniform'
        mskSPvol = sum(mskSP(:)); %volume of mask used to select start points
        NSP = round(gladpar.spPerc*mskSPvol/100);
        if ~exist(sprintf('%s_%s%d_spPerc%d_nsp%d.mat',tag,gladpar.spMsk_name,gladpar.spMsk_ind,gladpar.spPerc,NSP),'file')
            [spIND,spINDid] = datasample(mind,NSP,'Replace',false);
            %spIND = mind(spINDid);
            save(sprintf('%s_%s%d_spPerc%d_nsp%d.mat',tag,gladpar.spMsk_name,gladpar.spMsk_ind,gladpar.spPerc,NSP),'spIND');
        else
            load(sprintf('%s_%s%d_spPerc%d_nsp%d.mat',tag,gladpar.spMsk_name,gladpar.spMsk_ind,gladpar.spPerc,NSP));
        end
        [sy,sx,sz] = ind2sub(n,sort(spIND,'ascend'));
end
 
%% Lagrangian pathlines
%variables if nt > 1
mpar.n = n';
h1 = 1; h2 = 1; h3 = 1;
mpar.h1 = h1.*ones(n(1),1);
mpar.h2 = h2.*ones(n(2),1);
mpar.h3 = h3.*ones(n(3),1);
[mpar.Xc,mpar.Yc,mpar.Zc] = getCellCenteredGrid(mpar.h1,mpar.h2,mpar.h3);
mpar.Grad = getCellCenteredGradMatrix({'ccn' 'ccn' 'ccn'},mpar.h1,mpar.h2,mpar.h3);
Mdis = -mpar.sigma.*mpar.Grad'*mpar.Grad;
%initialize streamlines:
nsp = length(sx);
%convert from matlab grid to cell-centered grid:
s1 = (sy-0.5).*h1; %i/y-axis
s2 = (sx-0.5).*h2; %j/x-axis
s3 = (sz-0.5).*h3; %k/z-axis
sp_123 = [s1,s2,s3];


pcur = sp_123; %current point i.e. list of current location in each streamline that hasn't been terminated
npoints = length(pcur); %keep track of the # of streamlines that have not yet been terminated 

SL = cell(1,nsp);
RHO_SL = cell(1,nsp);
AUGSPD_SL = cell(1,nsp);
SPD_SL = cell(1,nsp);
DSPD_SL = cell(1,nsp);
T_SL = cell(1,nsp); %store timestep
FLX_SL = cell(1,nsp); %store flux
dR_SL = cell(1,nsp); %store drho/dt
PE_SL = cell(1,nsp); % store peclet

% Initialize masks
maskCluster = model_mask;
maskRho = model_mask;
maskSpeed = model_mask;
maskAugSpeed = model_mask;
maskTime = model_mask;
maskFlux = model_mask;
maskdRho = model_mask;
maskPe = model_mask;

% make sure don't use sl-info from previous timestep
File1 = fullfile(cd, 'pl_cur.mat');
File2 = fullfile(cd, 'pli_array.mat');
File3 = fullfile(cd, 'pl_centroid_array.mat');
if exist(File1,'file')
    delete(File1);
end
if exist(File2,'file')
    delete(File2);
end
if exist(File3,'file')
    delete(File3);
end

%initiate pathlines
pl = NaN(npoints,3,length(ti:tj:tf)*nt*gladpar.nEulStep+1);
plrho = NaN(npoints,length(ti:tj:tf)*nt*gladpar.nEulStep+1);
plspd = NaN(npoints,length(ti:tj:tf)*nt*gladpar.nEulStep+1);
pldspd = NaN(npoints,length(ti:tj:tf)*nt*gladpar.nEulStep+1);
plaugspd = NaN(npoints,length(ti:tj:tf)*nt*gladpar.nEulStep+1);
plt = NaN(npoints,length(ti:tj:tf)*nt*gladpar.nEulStep+1);
plflx = NaN(npoints,length(ti:tj:tf)*nt*gladpar.nEulStep+1);
pldr =  NaN(npoints,length(ti:tj:tf)*nt*gladpar.nEulStep+1);

step = 1;
xt = pcur; %current position in rcl orientation
xt = max([h1*0.5,h2*0.5,h3*0.5],min(xt,[h1*(n(1)-.5001),h2*(n(2)-.5001),h3*(n(3)-.5001)]));

%% running pathlines

for t1 = ti:tj:tf
    U = load(sprintf('%s/u0_%s_%d_%d_t_%d.mat',ppar.out_dir,tag,t1,t1+tj,(t1-ti)/tj + 1), 'u');
    U = reshape(U.u,[],nt);
    
    if strcmp(gladpar.RD,'R')
        if t1 == ti
            RHO = load(sprintf('%s/rho_%s_%d_t_0.mat',ppar.out_dir,tag,ti));
        else
            RHO = load(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',ppar.out_dir,tag,t1-tj,t1,(t1-ti)./tj));
        end
        RHO = RHO.rho_n;
        if nt > 1
            RHO_t = [RHO, advecDiff(RHO,U(:),nt,mpar.dt,mpar)];
        else
            RHO_t = RHO;
        end
    end
    
    for t2 = 1:nt
        TIND = ((t1 - ti)/tj)*nt + t2;
        T = t1+(t2-1)*(tj/nt);
        fprintf('t = %d (t1 = %d, t2 = %d -> T = %.3f)\n',TIND,t1,t2,T);
        
        switch gladpar.RD
            case 'D'
                d = getData(tag,round(T),'none');
                if dpar.smooth>0
                    d = affine_diffusion_3d(d,dpar.smooth,0.1,1,1);
                end
            case 'R'
                d = reshape(RHO_t(:,t2),n);
        end
        
        if gladpar.minIm0
            d = d-min(d(:));
        end
        
        %compute streamlines
        u = reshape(U(:,t2),[],3);
        v1 = reshape(u(:,1),n);
        v2 = reshape(u(:,2),n);
        v3 = reshape(u(:,3),n);
        
        %add eps to d so can take log(d) and not log(0)
        [w2,w1,w3] = gradient(log(d+2*eps));
        u1 = v1 - mpar.sigma.*w1;
        u2 = v2 - mpar.sigma.*w2;
        u3 = v3 - mpar.sigma.*w3;
        
        du = mpar.sigma*[w1(:),w2(:),w3(:)];
        
        switch gladpar.flw_type
            case 'vel'
                a1=u1;
                a2=u2;
                a3=u3;
            case 'flw'
                a1=u1.*d;
                a2=u2.*d;
                a3=u3.*d;
        end
           
        [Gx, Gy, Gz] = gradient(d);
        speed = reshape(sqrt(sum(u.^2,2)),n);
        dspeed = reshape(sqrt(sum(du.^2,2)),n);
        speedAug = sqrt(u1.^2 + u2.^2 + u3.^2);
        img_flow = speed.*d;
        drdt_dif = Mdis*d(:);
        drdt_ad1 = Gx.*reshape(u(:,2),n) + Gy.*reshape(u(:,1),n) + Gz.*reshape(u(:,3),n);
        DIVu = divergence(x,y,z,reshape(u(:,2),n),reshape(u(:,1),n),reshape(u(:,3),n));
        drdt_ad2 = d(:).*DIVu(:);
        drdt = reshape(drdt_dif - (drdt_ad1(:)+drdt_ad2),n);
        
        %update first step of density and speed
        if step == 1

            pl(:,:,1) = xt;
            plrho(:,1) = d(sub2ind(n,sy,sx,sz));
            plspd(:,1) = speed(sub2ind(n,sy,sx,sz));
            pldspd(:,1) = dspeed(sub2ind(n,sy,sx,sz));
            plaugspd(:,1) = speedAug(sub2ind(n,sy,sx,sz));
            plt(:,1) = step;%tstep;
            plflx(:,1) = img_flow(sub2ind(n,sy,sx,sz));
            pldr(:,1) = drdt(sub2ind(n,sy,sx,sz));
        end
        
        switch gladpar.cutoff_str
            case 'min'
                conf.conc = 1;
                conf.speed = 1;
            case 'max'
                conf.conc = mean(d(d>0)) + std(d(d>0));
                conf.speed = mean(speed(speed>0)) + std(speed(speed>0));              
            case 'mean'
                conf.conc = mean(d(d>0));
                conf.speed = mean(speed(speed>0));
        end
        
        %vector field to be integrated in order to compute streamlines
        V = [a1(:),a2(:),a3(:)];
        V(msk==0,:) = 0; %don't want to move outside of the masked ROI
        u(msk==0,:) = 0;
        
        for Estep = 1:gladpar.nEulStep
            step = step + 1;
            V_interp = interp_vel(V,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            u_interp = interp_vel(u,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            du_interp = interp_vel(du,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            
            conc_interp = interpF(d,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            spdAug_interp = sqrt(sum(V_interp.^2,2));
            spd_interp = sqrt(sum(u_interp.^2,2));
            dspd_interp = sqrt(sum(du_interp.^2,2));
            flx_interp = interpF(img_flow,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            dr_interp = interpF(drdt,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            
            conc_conf_interp = conc_interp./conf.conc;
            spd_conf_interp = spd_interp./conf.speed;
            
            thrsh_ind = find(conc_conf_interp > gladpar.thresholds.conc & spd_conf_interp > gladpar.thresholds.speed);

            if isempty(thrsh_ind)
                break
            else
                switch gladpar.XT
                    case'T'
                        xt(thrsh_ind,:) = xt(thrsh_ind,:) + gladpar.mdt.*V_interp(thrsh_ind,:);%so keep all current locations in case can take a step later when conc gets there
                    case 'X'
                        xt(thrsh_ind,:) = xt(thrsh_ind,:) + gladpar.mdt.*(V_interp(thrsh_ind,:)./spd_interp);%so keep all current locations in case can take a step later when conc gets there
                end
                %make sure it stays in bounds:
                xt = max([h1*0.5,h2*0.5,h3*0.5],min(xt,[h1*(n(1)-.5001),h2*(n(2)-.5001),h3*(n(3)-.5001)]));
                pl(thrsh_ind,:,step) = xt(thrsh_ind,:);
                plrho(thrsh_ind,step)  = conc_interp(thrsh_ind);
                plspd(thrsh_ind,step) = spd_interp(thrsh_ind);
                pldspd(thrsh_ind,step) = dspd_interp(thrsh_ind);
                plaugspd(thrsh_ind,step) = spdAug_interp(thrsh_ind);
                plt(thrsh_ind,step) =  step;
                plflx(thrsh_ind,step) = flx_interp(thrsh_ind);
                pldr(thrsh_ind,step) = dr_interp(thrsh_ind);
            end
        end
    end
end


%% select qualified pathlines
pli = 0; %counter for added pathlines
for sli = 1:npoints
    pl_cur = squeeze(pl(sli,:,:))';
    aind = any(~isnan(pl_cur),2); %1 if row is not NaN, 0 if row is NaN
    pl_cur = pl_cur(aind,:);
    %check that unique sline has more than 1 point (remove not-move coordinates)
    [pl_cur,ia,ic] = unique(pl_cur,'rows','stable');
    if size(pl_cur,1)>gladpar.pln
        pli = pli + 1;
        SL{pli} = pl_cur;
        
        plr_cur = plrho(sli,aind)';
        pls_cur = plspd(sli,aind)';
        plds_cur = pldspd(sli,aind)';
        plas_cur = plaugspd(sli,aind)';
        plt_cur = plt(sli,aind)';
        plflx_cur = plflx(sli,aind)';
        pldr_cur = pldr(sli,aind)';
        
        RHO_SL{pli} = plr_cur(ia);
        SPD_SL{pli} = pls_cur(ia);
        DSPD_SL{pli} = plds_cur(ia);
        AUGSPD_SL{pli} = plas_cur(ia);
        T_SL{pli} = plt_cur(ia);
        FLX_SL{pli} = plflx_cur(ia);
        dR_SL{pli} = pldr_cur(ia);
        PE_SL{pli} = pls_cur(ia)./(plds_cur(ia)+eps);
    end
end

%remove empty cell spaces from pathlines that were removed:
SL = SL(1:pli);
RHO_SL = RHO_SL(1:pli);
SPD_SL = SPD_SL(1:pli);
DSPD_SL = DSPD_SL(1:pli);
AUGSPD_SL = AUGSPD_SL(1:pli);
T_SL = T_SL(1:pli);
FLX_SL = FLX_SL(1:pli);
dR_SL = dR_SL(1:pli);
PE_SL = PE_SL(1:pli);

sl_euc = cellfun(@(x) sqrt(sum((x(end,:)-x(1,:)).^2)),SL);% getcell array with euclidean length of sl between first and last point
%figure, histogram(sl_euc),title('sl-euc'),axis tight;
SL2 = SL(sl_euc>gladpar.sl_tol);%only keep streamlines whose Euclidean length between first and last points is larger than the threshold
rstream2 = RHO_SL(sl_euc>gladpar.sl_tol);
sstream2 = SPD_SL(sl_euc>gladpar.sl_tol);
dsstream2 = DSPD_SL(sl_euc>gladpar.sl_tol);
asstream2 = AUGSPD_SL(sl_euc>gladpar.sl_tol);
tstream2 = T_SL(sl_euc>gladpar.sl_tol);
fstream2 = FLX_SL(sl_euc>gladpar.sl_tol);
dstream2 = dR_SL(sl_euc>gladpar.sl_tol);
pestream2 = PE_SL(sl_euc>gladpar.sl_tol);

fprintf(' # of start points = %d\n # of effective pathlines after pathline-number (pln) threshold = %d \n # of effective pathlines after Euclidean dist (sl_tol) threshold = %d\n',npoints,pli,length(SL2))
pl_cur = cellfun(@(x) x(:,[2,1,3]),SL2,'UniformOutput',false);


switch type
    case 's'
        save('pl_cur.mat','pl_cur');

        %save(sprintf('%s/%s/pl_cur.mat',ppar.out_dir,outdir),'pl_cur');
        %save(sprintf('%s/%s/pl.mat',ppar.out_dir,outdir),'pl');
        %save(sprintf('%s/%s/plspd.mat',ppar.out_dir,outdir),'plspd');
        %save(sprintf('%s/%s/plrho.mat',ppar.out_dir,outdir),'plrho');
        %save(sprintf('%s/%s/plaugspd.mat',ppar.out_dir,outdir),'plaugspd');
        %save(sprintf('%s/%s/plt.mat',ppar.out_dir,outdir),'plt');
        %save(sprintf('%s/%s/plflx.mat',ppar.out_dir,outdir),'plflx');
        %save(sprintf('%s/%s/pldr.mat',ppar.out_dir,outdir),'pldr');
        %% === run QB.py Here ===
        fprintf('\n======= waiting to run run_dipyQB_pl.py =======\n')
        fprintf('\nInstructions:\n\npl_cur.mat has been saved at the current directory.\nDirectly run run_dipyQB_pl.py with Python also in this directory.\nResults will be saved automatically.\nThen come back to Matlab and press any key to continue.\n\n')
        pause
        fprintf('...Matlab code sucessfully continues...\n')

        %% add start points of pathlines considered for clustering analysis
        load(fullfile(cd,'pli_array.mat'));
        load(fullfile(cd,'pl_centroid_array.mat'));
        % only want clusters with more than tol # of sls:
        clus_length = cellfun('size',pli_array,2)'; 
        fprintf(' # of original clusters = %d\n',length(clus_length))
        
        sli = pli_array(clus_length>gladpar.clus_tol);
        sl_centroid = pl_centroid_array(clus_length>gladpar.clus_tol,:,:);
        nclus = size(sli,2);%# of clusters that made the cutoff
        clus_length = cellfun('size',sli,2)'; 
        fprintf(' # of clusters after cluster-length (clus_tol) threshold = %d\n',length(clus_length))
        [~,clus_size_rank] = sort(clus_length);

        %only keep largest N clusters where N = clus_cutoff:
        if gladpar.clus_cutoff>0
            if nclus > gladpar.clus_cutoff
                ind_tmp = zeros(nclus,1);
                ind_tmp(clus_size_rank(end-gladpar.clus_cutoff+1:end)) = 1;
                %% to keep clusters in same order that they were returned from QB alg:
                sli = sli(ind_tmp==1);
                sl_centroid = sl_centroid(ind_tmp==1,:,:);
                nclus = size(sli,2);%# of clusters that made the cutoff
                fprintf('nclus=%d\n',nclus)
                clus_length = cellfun('size',sli,2)';
                [~,clus_size_rank] = sort(clus_length);
            end
        end
        fprintf(' # of clusters after max-cluster-number (clus_cutoff) threshold = %d\n',nclus)

        sl_centroid = double(sl_centroid);
        slc = cell(1,nclus);
        for ind = 1:nclus
            slc{1,ind}=squeeze(sl_centroid(ind,:,:));
        end
        [B,I] = sort(squeeze(sl_centroid(:,1,1)));%sl_centroid is nclus x npoints x 3

        % initialize temporary masks:
        c = zeros(n);%clusters
        r = zeros(n);%density
        s = zeros(n);%speed
        ds = zeros(n);%diffusive speed
        as = zeros(n);%augspeed
        tt = zeros(n);%time step
        ff = zeros(n);%flux
        drt = zeros(n);%dr/dt

        % make cluster masks
        %getting clustered pathline start points
        spTMP = zeros(n);
        for ind_clus = 1:nclus
            slines_tmp = pl_cur([sli{I(ind_clus)}]+1);
            rlines_tmp = rstream2([sli{I(ind_clus)}]+1);
            spdlines_tmp = sstream2([sli{I(ind_clus)}]+1);
            dspdlines_tmp = dsstream2([sli{I(ind_clus)}]+1);
            aspdlines_tmp = asstream2([sli{I(ind_clus)}]+1);
            tlines_tmp = tstream2([sli{I(ind_clus)}]+1);
            flines_tmp = fstream2([sli{I(ind_clus)}]+1);
            dlines_tmp = dstream2([sli{I(ind_clus)}]+1);

            %getting start points in world coords
            spx = cellfun(@(x) x(1,1),pl_cur([sli{I(ind_clus)}]+1));
            spy = cellfun(@(x) x(1,2),pl_cur([sli{I(ind_clus)}]+1));
            spz = cellfun(@(x) x(1,3),pl_cur([sli{I(ind_clus)}]+1));
            %getting start points in matlab coords
            spm1 = round(spy./h1 + 0.5);
            spm2 = round(spx./h2 + 0.5);
            spm3 = round(spz./h3 + 0.5);
            %save to spTMP
            indtmp = sub2ind(n,spm1,spm2,spm3);
            spTMP(indtmp) = 1;

            n_slines = size(slines_tmp,2);
            for ind_line = 1:n_slines
                %convert back to MATLAB grid
                sline = round((slines_tmp{ind_line}./[h1,h2,h3])+0.5);
                [sline,ia,~] = unique(sline, 'rows', 'stable');
                rsl = rlines_tmp{ind_line}(ia);
                ssl = spdlines_tmp{ind_line}(ia);
                dssl = dspdlines_tmp{ind_line}(ia);
                assl = aspdlines_tmp{ind_line}(ia);
                tsl = tlines_tmp{ind_line}(ia);
                fsl = flines_tmp{ind_line}(ia);
                dsl = dlines_tmp{ind_line}(ia);

                subs_1 = sline(:,2);
                subs_2 = sline(:,1);
                subs_3 = sline(:,3);
                subs_1(subs_1 < 1) = 1;
                subs_2(subs_2 < 1) = 1;
                subs_3(subs_3 < 1) = 1;
                subs_1(subs_1 > n(1)) = n(1);
                subs_2(subs_2 > n(2)) = n(2);
                subs_3(subs_3 > n(3)) = n(3);
                inds = sub2ind(n, subs_1, subs_2, subs_3);

                c(inds) = ind_clus*1;
                r(inds) = rsl;
                s(inds) = ssl;
                ds(inds) = dssl;
                as(inds) = assl;
                tt(inds) = tsl;
                ff(inds) = fsl;
                drt(inds) = dsl;
            end
            %centroid mask:
            if gladpar.do_centroid
                %centroid_tmp = round(slc{clus_size_rank(ind_clus)});
                centroid_tmp = round(slc{I(ind_clus)});
                centroid_tmp = unique(centroid_tmp,'rows', 'stable');
                subs_1 = centroid_tmp(:,2);
                subs_2 = centroid_tmp(:,1);
                subs_3 = centroid_tmp(:,3);
                subs_1(subs_1 < 1) = 1;
                subs_2(subs_2 < 1) = 1;
                subs_3(subs_3 < 1) = 1;
                subs_1(subs_1 > n(1)) = n(1);
                subs_2(subs_2 > n(2)) = n(2);
                subs_3(subs_3 > n(3)) = n(3);
                inds = sub2ind(n, subs_1, subs_2, subs_3);
                cent_tmp(inds) = ind_clus*1;
            end
        end

        if gladpar.do_masked
            %Added to ensure nothing outside of masked region
            c(~msk) = 0;
            r(~msk) = 0;
            s(~msk) = 0;
            ds(~msk) = 0;
            as(~msk) = 0;
            tt(~msk) = 0;
            ff(~msk) = 0;
            drt(~msk) = 0;
        end
        % only save within brain mask
        sp_ind = find(strcmp('brain',{ppar.sp_mask_opts(:).name}));
        mskROI = nii2mat(ppar.sp_mask_opts(sp_ind).path,dpar.x_range,dpar.y_range,dpar.z_range);
        msk_brain = mskROI>1;
        if dpar.do_resize
            msk_brain = resizeMatrix(double(msk_brain),round(dpar.size_factor.*size(msk_brain)),'linear');
            msk_brain(msk_brain~=1) = 0;
            c(msk_brain==0) = 0;
            r(msk_brain==0) = 0;
            s(msk_brain==0) = 0;
            ds(msk_brain==0) = 0;
            as(msk_brain==0) = 0;
            tt(msk_brain==0) = 0;
            ff(msk_brain==0) = 0;
            drt(msk_brain==0) = 0;
            
            Pe = s./ds;
            Pe(isnan(Pe)) = 0;
            Pe(Pe==Inf) = 0;
            
            maskCluster.img(round(dpar.x_range(1)*dpar.size_factor):round(dpar.x_range(end)*dpar.size_factor),round(dpar.y_range(1)*dpar.size_factor):round(dpar.y_range(end)*dpar.size_factor),round(dpar.z_range(1)*dpar.size_factor):round(dpar.z_range(end)*dpar.size_factor)) = c;
            maskRho.img(round(dpar.x_range(1)*dpar.size_factor):round(dpar.x_range(end)*dpar.size_factor),round(dpar.y_range(1)*dpar.size_factor):round(dpar.y_range(end)*dpar.size_factor),round(dpar.z_range(1)*dpar.size_factor):round(dpar.z_range(end)*dpar.size_factor)) = r;
            maskSpeed.img(round(dpar.x_range(1)*dpar.size_factor):round(dpar.x_range(end)*dpar.size_factor),round(dpar.y_range(1)*dpar.size_factor):round(dpar.y_range(end)*dpar.size_factor),round(dpar.z_range(1)*dpar.size_factor):round(dpar.z_range(end)*dpar.size_factor)) = s;
            maskAugSpeed.img(round(dpar.x_range(1)*dpar.size_factor):round(dpar.x_range(end)*dpar.size_factor),round(dpar.y_range(1)*dpar.size_factor):round(dpar.y_range(end)*dpar.size_factor),round(dpar.z_range(1)*dpar.size_factor):round(dpar.z_range(end)*dpar.size_factor)) = as;
            maskTime.img(round(dpar.x_range(1)*dpar.size_factor):round(dpar.x_range(end)*dpar.size_factor),round(dpar.y_range(1)*dpar.size_factor):round(dpar.y_range(end)*dpar.size_factor),round(dpar.z_range(1)*dpar.size_factor):round(dpar.z_range(end)*dpar.size_factor)) = tt;
            maskFlux.img(round(dpar.x_range(1)*dpar.size_factor):round(dpar.x_range(end)*dpar.size_factor),round(dpar.y_range(1)*dpar.size_factor):round(dpar.y_range(end)*dpar.size_factor),round(dpar.z_range(1)*dpar.size_factor):round(dpar.z_range(end)*dpar.size_factor)) = ff;
            maskdRho.img(round(dpar.x_range(1)*dpar.size_factor):round(dpar.x_range(end)*dpar.size_factor),round(dpar.y_range(1)*dpar.size_factor):round(dpar.y_range(end)*dpar.size_factor),round(dpar.z_range(1)*dpar.size_factor):round(dpar.z_range(end)*dpar.size_factor)) = drt;
            maskPe.img(round(dpar.x_range(1)*dpar.size_factor):round(dpar.x_range(end)*dpar.size_factor),round(dpar.y_range(1)*dpar.size_factor):round(dpar.y_range(end)*dpar.size_factor),round(dpar.z_range(1)*dpar.size_factor):round(dpar.z_range(end)*dpar.size_factor)) = Pe;
        else
            c(msk_brain==0) = 0;
            r(msk_brain==0) = 0;
            s(msk_brain==0) = 0;
            ds(msk_brain==0) = 0;
            as(msk_brain==0) = 0;
            tt(msk_brain==0) = 0;
            ff(msk_brain==0) = 0;
            drt(msk_brain==0) = 0;
            
            Pe = s./ds;
            Pe(isnan(Pe)) = 0;
            Pe(Pe==Inf) = 0;
            
            maskCluster.img(dpar.x_range,dpar.y_range,dpar.z_range) = c;
            maskRho.img(dpar.x_range,dpar.y_range,dpar.z_range) = r;
            maskSpeed.img(dpar.x_range,dpar.y_range,dpar.z_range) = s;
            maskAugSpeed.img(dpar.x_range,dpar.y_range,dpar.z_range) = as;
            maskTime.img(dpar.x_range,dpar.y_range,dpar.z_range) = tt;
            maskFlux.img(dpar.x_range,dpar.y_range,dpar.z_range) = ff;
            maskdRho.img(dpar.x_range,dpar.y_range,dpar.z_range) = drt;
            maskPe.img(dpar.x_range,dpar.y_range,dpar.z_range) = Pe;
        end
        
        
        % save to nifty (view in Amira later)
        
        %save_untouch_nii(maskCluster,sprintf('%s/%s/%s_LagClusters_E%02d_%02d_%s_%s.nii',ppar.out_dir,outdir,tag,ti,tf+tj,paper_fig_str,date_str));
        %save_untouch_nii(maskRho,sprintf('%s/%s/%s_LagRho_E%02d_%02d_%s_%s.nii',ppar.out_dir,outdir,tag,ti,tf+tj,paper_fig_str,date_str));
        save_untouch_nii(maskSpeed,sprintf('%s/%s/%s_LagSpeed_E%02d_%02d_%s_%s.nii',ppar.out_dir,outdir,tag,ti,tf+tj,paper_fig_str,date_str));
        %save_untouch_nii(maskAugSpeed,sprintf('%s/%s/%s_LagAugSpeed_E%02d_%02d_%s_%s.nii',ppar.out_dir,outdir,tag,ti,tf+tj,paper_fig_str,date_str));
        %save_untouch_nii(maskTime,sprintf('%s/%s/%s_LagTime_E%02d_%02d_%s_%s.nii',ppar.out_dir,outdir,tag,ti,tf+tj,paper_fig_str,date_str));
        %save_untouch_nii(maskFlux,sprintf('%s/%s/%s_LagRhoFlux_E%02d_%02d_%s_%s.nii',ppar.out_dir,outdir,tag,ti,tf+tj,paper_fig_str,date_str));
        %save_untouch_nii(maskdRho,sprintf('%s/%s/%s_LagdRhodt_E%02d_%02d_%s_%s.nii',ppar.out_dir,outdir,tag,ti,tf+tj,paper_fig_str,date_str));
        save_untouch_nii(maskPe,sprintf('%s/%s/%s_LagPe_E%02d_%02d_%s_%s.nii',ppar.out_dir,outdir,tag,ti,tf+tj,paper_fig_str,date_str));
        
        %fprintf('For Speed Map, the next step is to visualize the speed nifty file in Amira.\nHowever, we can still plot in Matlab\n')
        fprintf('Speed and Peclet Map in nifty format saved in %s/%s\n\n',ppar.out_dir,outdir)
        
        %% Visualization
        x = 1:n(1);
        y = 1:n(2);
        z = 1:n(3);
        figure,
        hs=slice(y,x,z,s,x,y,z); 
        set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
        alpha('color'),alphamap(linspace(0,1,100))
        title(sprintf('Test: tag = %s, Speed Map',tag),'Fontsize',20)
        grid off, box off, axis image
        xlabel('x-axis')
        ylabel('y-axis')
        zlabel('z-axis')
        colormap(jet)
        caxis([0,0.5])
        view([-188.3500   13.7463])
        set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')
        colorbar, grid on,

        saveas(gcf, sprintf('%s/%s/%s_LagSpeed_E%02d_%02d.png',ppar.out_dir,outdir,tag,ti,tf+tj)); 
        
        figure,
        hs=slice(y,x,z,Pe,x,y,z); 
        set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
        alpha('color'),alphamap(linspace(0,1,100))
        title(sprintf('Test: tag = %s, Peclet Map',tag),'Fontsize',20)
        grid off, box off, axis image
        xlabel('x-axis')
        ylabel('y-axis')
        zlabel('z-axis')
        colormap(jet)
        caxis([0,200])
        view([-188.3500   13.7463])
        set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')
        colorbar, grid on,

        saveas(gcf, sprintf('%s/%s/%s_LagPe_E%02d_%02d.png',ppar.out_dir,outdir,tag,ti,tf+tj)); 
    case 'v'

        outversion = sprintf('%s_%s',paper_fig_str,date_str);
        outdir = sprintf('LPPA_%s',outversion);

        if ~exist(sprintf('%s/%s',ppar.out_dir,outdir),'dir')
            fprintf('%s/%s\n Directory does not exist :(\n',ppar.out_dir,outdir);
            return
        else
            fprintf('%s: %s Directory exists :)\n',tag,outdir);
        end

        %sstream = sstream2;
        fprintf('Total original %d pathlines\n',length(SL2)); 

        %convert from cell-centered grid to matlab grid:
        SL = cellfun(@(x) [x(:,1)./h1+0.5,x(:,2)./h2+0.5,x(:,3)./h3+0.5],SL2,'UniformOutput',false); 
        clear SL2;
        
        dispnor = cellfun(@(x) (x(end,:)-x(1,:))/norm(x(end,:)-x(1,:)),SL,'UniformOutput',false); 
        disp = cellfun(@(x) x(end,:)-x(1,:),SL,'UniformOutput',false); 
        startp = cellfun(@(x) x(1,:),SL,'UniformOutput',false); 
        endp = cellfun(@(x) x(end,:),SL,'UniformOutput',false); 
        
        PATH.NPtsInPath = cellfun(@(x) size(x,1),SL); % number of points in each pathline
        PATH.LengthInPath = cellfun(@(x) sum(sqrt(sum(diff(x).^2,2))),SL); % total length of path in each pathline
        PATH.disp = reshape([disp{:}]',3,[])'; % displacement field
        PATH.dispnor = reshape([dispnor{:}]',3,[])'; % normalized displacement field
        PATH.startp = reshape([startp{:}]',3,[])'; % startp points
        PATH.endp = reshape([endp{:}]',3,[])'; % endp points
        PATH.displen = sqrt(PATH.disp(:,1).^2+PATH.disp(:,2).^2+PATH.disp(:,3).^2); % displacement length in each pathline

        anato = load_untouch_nii(ppar.anato); 
        anato = anato.img(dpar.x_range,dpar.y_range,dpar.z_range);
        sp_ind = find(strcmp('brain',{ppar.sp_mask_opts(:).name}));
        mskROI = nii2mat(ppar.sp_mask_opts(sp_ind).path,dpar.x_range,dpar.y_range,dpar.z_range);
        msk_brain = mskROI>1;
        
        if dpar.do_resize
           anato = resizeMatrix(anato,round(dpar.size_factor.*size(anato)),'linear');
           msk_brain = resizeMatrix(double(msk_brain),round(dpar.size_factor.*size(msk_brain)),'linear');
           msk_brain(msk_brain~=1) = 0;
        end
        anato(~msk) = 0;
        anato(msk_brain==0) = 0;
        indb = find(msk_brain(sub2ind(n,PATH.startp(:,1),round(PATH.startp(:,2)),PATH.startp(:,3)))==1); % index of those starting within brain
        
        % save to vtk
        SL_tmp = SL(indb);
        PE_tmp = pestream2(indb);
        S_tmp = sstream2(indb);
        vtkwrite_pathlines(sprintf('%s/%s/%s_pathlines_lentol_%.1f_jp_%d_%s.vtk',ppar.out_dir,outdir,tag,gladpar.sl_tol,gladpar.jp,outversion),'polydata','lines',SL_tmp(1:gladpar.jp:end));
        vtkwrite_spdlines(sprintf('%s/%s/%s_Pelines_lentol_%.1f_jp_%d_%s.vtk',ppar.out_dir,outdir,tag,gladpar.sl_tol,gladpar.jp,outversion),'polydata','lines',SL_tmp(1:gladpar.jp:end),PE_tmp(1:gladpar.jp:end));
        vtkwrite_spdlines(sprintf('%s/%s/%s_Spdlines_lentol_%.1f_jp_%d_%s.vtk',ppar.out_dir,outdir,tag,gladpar.sl_tol,gladpar.jp,outversion),'polydata','lines',SL_tmp(1:gladpar.jp:end),S_tmp(1:gladpar.jp:end));
        vtkwrite(sprintf('%s/%s/%s_disp_lentol_%.2f_%s.vtk',ppar.out_dir,outdir,tag,gladpar.sl_tol,outversion), 'structured_grid', PATH.startp(indb,1), PATH.startp(indb,2), PATH.startp(indb,3), ... 
            'vectors', 'vector_field', PATH.disp(indb,1), PATH.disp(indb,2), PATH.disp(indb,3));
        vtkwrite(sprintf('%s/%s/%s_anato_%s.vtk',ppar.out_dir,outdir,tag,outversion), 'structured_points', 'mask', anato);
        
        fprintf('Pathlines and Flux vectors in vtk format saved in %s/%s\n\n',ppar.out_dir,outdir)
        %% NCA
        %{
        [xx, yy, zz] = meshgrid(1:n(2),1:n(1),1:n(3));

        startpind = sub2ind(n,round(PATH.startp(:,1)),round(PATH.startp(:,2)),round(PATH.startp(:,3)));

        Mean = zeros(length(PATH.displen),1); % mean in the neighbor
        WMean = Mean; % weighted mean in the neighbor
        STD = Mean; % L_2^2 distance to 1
        Np = Mean; % number of all points in neighbored paths
        Avepathlen = Mean; % mean path length in neighborhood
        %AveMaxpathspd = Mean; % average max speed  
        NNnum = zeros(length(PATH.displen),1); % number of neighbors

        for sp = 1:length(PATH.displen)
            ctr = round(PATH.startp(sp,:));
            binaryMap = (yy-ctr(1)).^2+(xx-ctr(2)).^2+(zz-ctr(3)).^2 <=gladpar.radius^2;
            indd = find(binaryMap); % startpoint ind in n1*n2*n3 coord.
            % pick neighboring points
            [~,Locb] = ismember(indd,startpind);
            Lia = Locb(Locb~=0);

            NNnum(sp) = length(Lia);
            costheta = dot(repmat(PATH.disp(sp,:),length(Lia),1),PATH.disp(Lia,:),2)/norm(PATH.disp(sp,:))./vecnorm(PATH.disp(Lia,:),2,2);
            WMean(sp) = sum(costheta.*PATH.displen(Lia))/sum(PATH.displen(Lia));
            Mean(sp) = mean(costheta);
            STD(sp) = sum((costheta-1).^2)/length(costheta);%std(costheta);
            Np(sp) = sum(PATH.NPtsInPath(Lia));
            Avepathlen(sp) = mean(PATH.LengthInPath(Lia));
            %AveMaxpathspd(sp) = mean(cellfun(@(x) max(x),sstream(Lia)));  

        end

        %%
        INDD = find(NNnum > gladpar.NNnum_tol & STD<gladpar.stdcut & Np>gladpar.Npcut & WMean>=gladpar.meancut(1) & Avepathlen>gladpar.Avepathlcut);
        INDD2 = find(NNnum > gladpar.NNnum_tol);

        %tmpind = INDD;
        %tmpind2 = setdiff(INDD2,INDD);
        %PATH.ADVind = tmpind;
        %PATH.DIFFind = tmpind2;
        %PATH.threshstr = gladpar.threshstr;

        %% further choose potential adv vectors
        maskadv = zeros(n);
        maskadv(sub2ind(n,round(PATH.startp(INDD,1)),round(PATH.startp(INDD,2)),round(PATH.startp(INDD,3))))=1;

        [xr,yr,zr] = meshgrid(-gladpar.maskdilate:gladpar.maskdilate,-gladpar.maskdilate:gladpar.maskdilate,-gladpar.maskdilate:gladpar.maskdilate);
        strel = (xr/gladpar.maskdilate).^2 + (yr/gladpar.maskdilate).^2 + (zr/gladpar.maskdilate).^2 <= 1;
        maskadvdia = imdilate(maskadv,strel);

        % find potential start points
        indp = find(maskadvdia(sub2ind(n,PATH.startp(:,1),round(PATH.startp(:,2)),round(PATH.startp(:,3))))==1 ...
            & maskadv(sub2ind(n,PATH.startp(:,1),round(PATH.startp(:,2)),round(PATH.startp(:,3))))==0);
        indpADV = find(maskadv(sub2ind(n,PATH.startp(:,1),round(PATH.startp(:,2)),round(PATH.startp(:,3))))==1);

        % set parameters
        %%
        Mean2 = -1*ones(length(PATH.displen),1); % mean in the neighbor
        WMean2 = Mean2; % weighted mean in the neighbor
        STD2 = Mean2; % L_2^2 distance to 1
        Np2 = Mean2; % number of all points in neighbored paths
        Avepathlen2 = Mean2; % mean path length in neighborhood
        %AveMaxpathspd2 = Mean2; % average max speed 
        pathl2 = Mean2; % path length in sp
        NNnum2 = -1*ones(length(PATH.displen),1); % number of neighbors

        for sp = 1:length(PATH.displen)
            if ~ismember(sp,indp)%sp is not in indp
                continue
            end
            ctr = round(PATH.startp(sp,:));
            binaryMap = (yy-ctr(1)).^2+(xx-ctr(2)).^2+(zz-ctr(3)).^2 <=gladpar.radius2^2;
            indd = find(binaryMap==1); % startpoint ind in n1*n2*n3 coord. 
            [~,Locb] = ismember(indd,indpADV); % only find neighbor in maskadv
            Lia = Locb(Locb~=0);

            NNnum2(sp) = length(Lia);
            costheta = dot(repmat(PATH.disp(sp,:),length(Lia),1),PATH.disp(Lia,:),2)/norm(PATH.disp(sp,:))./vecnorm(PATH.disp(Lia,:),2,2);
            %costheta(costheta<=-0.75) = 1; % remove those are too opposite trend
            WMean2(sp) = sum(costheta.*PATH.displen(Lia))/sum(PATH.displen(Lia));
            Mean2(sp) = mean(costheta);
            STD2(sp) = sum((costheta-1).^2)/length(costheta);%std(costheta);
            Np2(sp) = sum(PATH.NPtsInPath(Lia));
            Avepathlen2(sp) = mean(PATH.LengthInPath(Lia));
            pathl2(sp) = PATH.LengthInPath(sp);
            %AveMaxpathspd2(sp) = mean(cellfun(@(x) max(x),sstream(Lia)));
        end

        %%
        INDD_dia = find(NNnum2>gladpar.NNnum_tol2 & STD2<gladpar.stdcut2 & Np2>gladpar.Npcut2 & WMean2>=gladpar.meancut2(1) & Avepathlen2>gladpar.Avepathlcut2 | pathl2>gladpar.pathlcut2);
        fprintf('After further dilate to add more ADV, %d vectors are added among %d candidates to %d already ADV vectors\n',length(INDD_dia),length(indp),length(INDD))

        tmpind = [INDD;INDD_dia];
        tmpind2 = setdiff(INDD2,tmpind);
        PATH.ADVind = tmpind;
        PATH.DIFFind = tmpind2;
        PATH.threshstr = sprintf('Overall separation: \n%s, \nandfurther dilate to select more ADV: \n%s',gladpar.threshstr,gladpar.threshstr2);
        
        [InDadv,~,~] = intersect(PATH.ADVind,indb);
        [InDdiff,~,~] = intersect(PATH.DIFFind,indb);
        % save
        vtkwrite(sprintf('%s/%s/%s_ADVdisp_%s.vtk',ppar.out_dir,outdir,tag,outversion), 'structured_grid', PATH.startp(InDadv,1), PATH.startp(InDadv,2), PATH.startp(InDadv,3), ... 
            'vectors', 'vector_field', PATH.disp(InDadv,1), PATH.disp(InDadv,2), PATH.disp(InDadv,3));

        vtkwrite(sprintf('%s/%s/%s_DIFFdisp_%s.vtk',ppar.out_dir,outdir,tag,outversion), 'structured_grid', PATH.startp(InDdiff,1), PATH.startp(InDdiff,2), PATH.startp(InDdiff,3), ... 
            'vectors', 'vector_field', PATH.disp(InDdiff,1), PATH.disp(InDdiff,2), PATH.disp(InDdiff,3));
        
        fprintf('Flux vectors in vtk format saved in %s/%s\n\n',ppar.out_dir,outdir)
        %}
end

end
























