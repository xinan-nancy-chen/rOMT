function [dr0Tx, rN] = get_drNdr0T(rho0,u,nt,dt,par,x)
%% Sensitivities of rho_{n} w.r.t rho_{0} transpose times vector
%  x - vector size (prod(n),1)

n                  = par.n;
%rho                = zeros(prod(n),nt+1);
Mdis               = -par.sigma*par.Grad'*par.Grad;
I                  = speye(prod(n));
%rho(:,1)           = rho0;
dr0Tx              = speye(prod(n))'*x;
u                  = reshape(u,3*prod(n),[]);

for i = nt+1:-1:2
    U1 = reshape(u(1:prod(n),i-1),n');
    U2 = reshape(u(prod(n)+1:2*prod(n),i-1),n');
    U3 = reshape(u(2*prod(n)+1:end,i-1),n');
    S  = dTrilinears(rho0,par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3, ...
                     par.h1(1),par.h2(1),par.h3(1));
    
    %rho(:,i-1)  = S'*((I - dt*Mdis)\rho(:,i));
    %rho(:,i-1)  = S'*pcg((I - dt*Mdis),rho(:,i));
    %% Sensitivities
    %dr0Tx       = S'*((I - dt*Mdis)\dr0Tx);   
    %dr0Tx       = S'*pcg((I - dt*Mdis),dr0Tx);
    [dr0TxI,pcgflag]       = pcg((I - dt.*Mdis),dr0Tx);
    if pcgflag ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNdroT.m >>> while finding drho_%dT/drho_0, pcg exit flag = %d',i,pcgflag)
    end
    dr0Tx = S'*dr0TxI;
end;

rN = rho0;%rN = rho(:,end);
