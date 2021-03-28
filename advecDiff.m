function [rho] = advecDiff(rho0,u,nt,dt, par)
%% rho(:,1...T) = transport (rho0 - initial mass density, 
%                                 u - velocity vector for nt time steps,
%                                     size(prod(n),nt)
%                                 nt - time steps)
%                   rho_N - density at the last time step

% Advection: Semi- Lagrangian  rho^(i,ad) = S(U^(i)*rho^(i-1);

% Dispersion: Implicitly: (I - dt*Div*D*Grad)*rho^(i) = rho^(i,ad)

% rho(:,0)  -----------> rho(:,1)  --------------> rho(:,2)
%             u(:,1)                   u(:,2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n                  = par.n;
rho                = zeros(prod(n),nt+1);
u                  = reshape(u,3*prod(n),nt);
Mdis               = - par.sigma*par.Grad'*par.Grad;  % change to add variable sigma
I                  = speye(prod(n));
rho(:,1)           = rho0;

for i = 2:nt+1
    U1 = reshape(u(1:prod(n),i-1),n');
    U2 = reshape(u(prod(n)+1:2*prod(n),i-1),n');
    U3 = reshape(u(2*prod(n)+1:end,i-1),n');
    
    S  = dTrilinears(rho(:,i-1),par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3,...
                     par.h1(1),par.h2(1),par.h3(1));
    
    %T=rho(:,i-1);x=par.Xc + dt*U1;y= par.Yc + dt*U2;z= par.Zc + dt*U3;
    % hx=1;hy=1;hz=1;
    
    % advection step 
    %rho(:,i)  = S*rho(:,i-1);
    % advection step with source:
    if par.add_source
        rho(:,i)  = S*(rho(:,i-1) + dt.*par.qex(:));
    else
        rho(:,i)  = S*rho(:,i-1);
    end
    
    % diffusion step
       
    [rho(:,i),pcgflag]  = pcg(I - dt*Mdis,rho(:,i));
    if pcgflag ~= 0
        warning('MATLAB:pcgExitFlag','Warning: advecDiff.m >>> while finding rho(:,%d), pcg exit flag = %d',i,pcgflag)
    end
    %rho(:,i)  = (I - dt*Mdis)\rho(:,i);
    
    %montageArray(reshape(rho(:,i),n')); 
    %pause(0.1);
end;

%rho  = squeeze(rho(:,2:end));
rho = rho(:,2:end);

