% This is the driver for the rOMT algorithm
% parameters of rOMT is defined in getPamams.m

addpath(genpath('utilities', 'Sensitivities','Inverse'));

tag = 'C294';

%%
% get parameters and set directories
[ppar,dpar,mpar] = getParams(tag);

if dpar.do_time_interp
  time_interp = 'gauss';
 else
   time_interp = 'none';
end

reInitializeU = 1; %1 if reinitialize u to 0 before each time step; 0 if not, unless first time step

if ~exist(ppar.out_dir,'dir')
    mkdir(ppar.out_dir)
end

fname=sprintf('%s/record.txt',ppar.out_dir);
if ~exist(sprintf('%s/record.txt',ppar.out_dir),'file')
    csvwrite_with_headers(fname,[0 0 0 0 0 0 0 0 0],{'time-ind','ti','tf','phi','mk','Ru','phiN','max(u)','toc'});
end

if dpar.mask_number~=0
    mask = nii2mat(ppar.data_mask_path,dpar.x_range,dpar.y_range,dpar.z_range);
    msk = zeros(size(mask));
    msk(mask>0) = 1;
else
    msk = zeros(dpar.true_size);
    msk(dpar.x_range-dpar.x_range(1)+1,dpar.y_range-dpar.y_range(1)+1,dpar.z_range-dpar.z_range(1)+1) = 1;   
end

if dpar.do_resize
   msk = resizeMatrix(msk,round(dpar.size_factor.*size(msk)),'linear');
   msk(msk~=1) = 0;
end

if dpar.dilate>0
    [xr,yr,zr] = meshgrid(-dpar.dilate:dpar.dilate,-dpar.dilate:dpar.dilate,-dpar.dilate:dpar.dilate);
    strel = (xr/dpar.dilate).^2 + (yr/dpar.dilate).^2 + (zr/dpar.dilate).^2 <= 1;
    msk = imdilate(msk,strel);
end

rho_n = getData(tag,dpar.first_time,time_interp);
if dpar.smooth>0
    rho_n = affine_diffusion_3d(rho_n,dpar.smooth,0.1,1,1);
end
m=min(rho_n(:));

rho_n(~msk) = m;
rho_n = rho_n(:);

if ~exist(sprintf('%s/rho_%s_%d_t_0.mat',ppar.out_dir,ppar.data_tag,dpar.first_time),'file')
    save(sprintf('%s/rho_%s_%d_t_0.mat',ppar.out_dir,ppar.data_tag,dpar.first_time),'rho_n');
end
fprintf('\n =============== rOMT Starts ===============\n')
fprintf('______________________________________________\n\n')
fprintf(' tag:\t\t%s\n dataset:\t%s\n sigma:\t\t%.4f\n gamma:\t\t%.4f\n beta:\t\t%.4f\n nt:\t\t%d\n dt:\t\t%.2f\n pcg:\t\t%d\n',ppar.data_tag,dpar.dataset_name,mpar.sigma,mpar.gamma,mpar.beta,mpar.nt,mpar.dt,mpar.niter_pcg)
fprintf(' mask:\t\t%d\n do_resize:\t%d\n resize_factor:\t%.2f\n start frame:\t%d\n end frame:\t%d\n frame jump:\t%d\n\n\n',dpar.mask_number,dpar.do_resize,dpar.size_factor,dpar.first_time,dpar.last_time+dpar.time_jump,dpar.time_jump)
%%
clear T; T = 0;
for tind = 1:length(dpar.first_time:dpar.time_jump:dpar.last_time)
    fprintf('tind = %d\n',tind)
    tic
   
    if exist(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',ppar.out_dir,ppar.data_tag,dpar.first_time+(tind-1)*dpar.time_jump,dpar.first_time+(tind-1)*dpar.time_jump+dpar.time_jump,tind),'file')==2 && exist(sprintf('%s/u0_%s_%d_%d_t_%d.mat',ppar.out_dir,ppar.data_tag,dpar.first_time+(tind-1)*dpar.time_jump,dpar.first_time+(tind-1)*dpar.time_jump+dpar.time_jump,tind),'file') == 2        
        rho_n = load(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',ppar.out_dir,ppar.data_tag,dpar.first_time+(tind-1)*dpar.time_jump,dpar.first_time+(tind-1)*dpar.time_jump+dpar.time_jump,tind));
        rho_n = rho_n.rho_n;
        u = load(sprintf('%s/u0_%s_%d_%d_t_%d.mat',ppar.out_dir,ppar.data_tag,dpar.first_time+(tind-1)*dpar.time_jump,dpar.first_time+(tind-1)*dpar.time_jump+dpar.time_jump,tind));
        u = reshape(u.u,[],mpar.nt);
        continue
    end

    if dpar.reinitR
        rho_0 = getData(tag,dpar.first_time+(tind-1)*dpar.time_jump,time_interp);
        if dpar.smooth>0
	    rho_0 = affine_diffusion_3d(rho_0,dpar.smooth,0.1,1,1);
        end
        rho_0(~msk) = m;
        rho_0 = rho_0(:);
    else
        rho_0 = rho_n(:);% rho_I(:);
    end

    %true final density
    rho_N = getData(tag,dpar.first_time+tind*dpar.time_jump,time_interp);
    if dpar.smooth>0
        rho_N = affine_diffusion_3d(rho_N,dpar.smooth,0.1,1,1);
    end
    rho_N(~msk) = m;
    
    par = paramInitFunc(dpar.true_size',mpar.nt,mpar.dt,mpar.sigma,mpar.add_source,mpar.gamma,mpar.beta,mpar.niter_pcg);
    par.drhoN     = rho_N(:);
    
    % initial guess for u:
    if tind == 1 || reInitializeU
    u = zeros(3*prod(par.n),par.nt);
    end
    
    %% Descent for u
    fprintf('\n =============== Descent on u ===============\n')
    fprintf('______________________________________________\n\n')
    fprintf('i.lsiter\tphi    \t      descent output\n')
    fprintf('________    ___________     __________________\n')
    [u,phi,dphi] = GNblock_u(rho_0,u,par.nt,par.dt,par);
    
    
    [phi,mk,phiN,Rho,Ru]  = get_phi(rho_0,u,par.nt,par.dt,par);
    rho_n = Rho(:,end);
    btoc = toc;
    T = T + btoc;
    
    dlmwrite(fname,[tind,dpar.first_time+(tind-1)*dpar.time_jump,dpar.first_time+(tind-1)*dpar.time_jump+dpar.time_jump,phi,mk,Ru,phiN,max(u(:)),btoc],'-append');
    
    save(sprintf('%s/u0_%s_%d_%d_t_%d.mat',ppar.out_dir,ppar.data_tag,dpar.first_time+(tind-1)*dpar.time_jump,dpar.first_time+(tind-1)*dpar.time_jump+dpar.time_jump,tind),'u');
    save(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',ppar.out_dir,ppar.data_tag,dpar.first_time+(tind-1)*dpar.time_jump,dpar.first_time+(tind-1)*dpar.time_jump+dpar.time_jump,tind),'rho_n');
    
    fprintf('tind = %d, max(u) = %5.4f\n',tind,max(u));
end
fprintf('\n =============== rOMT Ends ===============\n')
fprintf('\n Elapsed Time: %s\n',datestr(seconds(T),'HH:MM:SS'))

