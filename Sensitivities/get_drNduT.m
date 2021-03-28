function [drNduT] = get_drNduT(rho_0,u,nt,dt,par,y)
%% Sensitivity of rho(:,end) w.r.t 'u' transpose times a vector x
%  y... vector length(prod(n),1)

% rho(:,0)  -----------> rho(:,1)  --------------> rho(:,2)
%             u(:,1)                   u(:,2)
%               S1                      S2

n                  = par.n;
u                  = reshape(u,3*prod(n),nt);
rho                = advecDiff(rho_0,u,nt,dt, par);
Mdis               = - par.sigma*par.Grad'*par.Grad;
I                  = speye(prod(n));
B                  = (I - dt*Mdis);

sensTx             = zeros(3*prod(n),nt);
%sens               = repmat(y,1,nt);
sens               = y;
Rho                = zeros(prod(n),nt+1);
Rho(:,1)           = rho_0; Rho(:,2:end) = rho;

for i = nt:-1:1
    U1 = reshape(u(1:prod(n),i),n');
    U2 = reshape(u(prod(n)+1:2*prod(n),i),n');
    U3 = reshape(u(2*prod(n)+1:end,i),n');
    
    [S,Tx,Ty,Tz]  = dTrilinears(Rho(:,i),par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3,...
                     par.h1(1),par.h2(1),par.h3(1));
    
    % Sensitivity:          
    %sensTx(:,i)              = dt*[Tx,Ty,Tz]'*(B'\sens(:,i));
    %sensTx(:,i)              = dt*[Tx,Ty,Tz]'*(pcg(B',sens(:,i)));
    %[sensI,pcgflag1]   = pcg(B',sens(:,i));
    [sensI,pcgflag1]   = pcg(B',sens);
    if pcgflag1 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNduT.m >>> while finding drho_%dT/du, pcg exit flag1 = %d',i,pcgflag1)
    end
    sensTx(:,i)              = dt*[Tx,Ty,Tz]'*sensI;
    
    if  i>1 
    %sens(:,1:i-1) = S'*(B'\sens(:,1:i-1));
    %sens(:,1:i-1) = S'*(pcg(B',:,1:i-1)));
    %[sensi,pcgflag2] = pcg(B',sens(:,1:i-1));
    [sensi,pcgflag2] = pcg(B',sens);
    if pcgflag2 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNduT.m >>> while finding drho_%dT/du, pcg exit flag2 = %d',i,pcgflag2)
    end
    %sens(:,1:i-1) = S'*sensi;
    sens = S'*sensi;
    end;
end;

drNduT = sensTx(:);