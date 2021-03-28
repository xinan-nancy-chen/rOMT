function [dr0x, rN] = get_drNdr0(rho0,u,nt,dt,par,x)
%% Sensitivities of rho_{n} w.r.t rho_{0}
%  x - vector size (prod(n),1)

n                  = par.n;
rho                = zeros(prod(n),nt+1);
Mdis               = -par.sigma*par.Grad'*par.Grad;
I                  = speye(prod(n));
rho(:,1)           = rho0;
dr0x               = x;
u                  = reshape(u,3*prod(n),[]);

for i = 2:nt+1
    U1 = reshape(u(1:prod(n),i-1),n');
    U2 = reshape(u(prod(n)+1:2*prod(n),i-1),n');
    U3 = reshape(u(2*prod(n)+1:end,i-1),n');
    S  = dTrilinears(rho(:,i-1),par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3,... 
                     par.h1(1),par.h2(1),par.h3(1));
    
    %rho(:,i)  = (I - dt*Mdis)\(S*rho(:,i-1));
    
    %% Sensitivities
    %dr0x     = (I - dt*Mdis)\(S*dr0x);
    [dr0x,pcgflag]     = pcg((I - dt*Mdis),(S*dr0x));
     if pcgflag ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNdr0.m >>> while finding drho_%d/drho_0, pcg exit flag = %d',i-1,pcgflag)
     end
    
end;

rN = rho(:,end);