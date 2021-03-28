function [drNdu] = get_drNdu(rho_0,u,nt,dt,par,x)
%% Sensitivity of rho(:,end) w.r.t 'u', full vector
% x... vector length(u,2), or 3*prod(n)*nt

% rho(:,0)  -----------> rho(:,1)  --------------> rho(:,2)
%             u(:,1)                   u(:,2)

n                  = par.n;
u                  = reshape(u,3*prod(n),nt);
rho                = zeros(prod(n),nt+1);
rho(:,1)           = rho_0;
Mdis               = -par.sigma*par.Grad'*par.Grad;
I                  = speye(prod(n));
B                  = (I - dt*Mdis);

X                  = reshape(x,prod(n)*3,nt);
sensx              = zeros(prod(n),nt);


for i = 1:nt
    U1 = reshape(u(1:prod(n),i),n');
    U2 = reshape(u(prod(n)+1:2*prod(n),i),n');
    U3 = reshape(u(2*prod(n)+1:end,i),n');
    
    [S,Tx,Ty,Tz]  = dTrilinears(rho(:,i),par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3,...
        par.h1(1),par.h2(1),par.h3(1));
    
    
    rho(:,i+1)  = S*rho(:,i);
    [rho(:,i+1),pcgflag1]  = pcg((I - dt*Mdis),rho(:,i+1));
    if pcgflag1 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNdu.m >>> while finding drho(:,%d)/du, pcg exit pcgflag1 = %d',i,pcgflag1)
    end
    %rho(:,i+1)  = (I - dt*Mdis)\rho(:,i+1);
    
    % Sensitivity:
    if  i>1
        for j = 1:i-1
            [sensx(:,j),pcgflag2] = pcg(B,(S*sensx(:,j)));
            %[sensx(:,1:i-1),pcgflag2] = pcg(B,(S*sensx(:,1:i-1)));
            %sensx(:,1:i-1) = B\(S*sensx(:,1:i-1));
            if pcgflag2 ~= 0
                warning('MATLAB:pcgExitFlag','Warning: get_drNdu.m >>> while finding drho(:,n)/du_%d, pcg exit pcgflag2 = %d',j,pcgflag2)
            end
        end
    end;
    [sensx(:,i),pcgflag3] = pcg(B,(dt*[Tx,Ty,Tz]*X(:,i)));
    if pcgflag3 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNdu.m >>> while finding drho(:,%d)/du, pcg exit pcgflag3 = %d',i,pcgflag3)
    end
    %sensx(:,i) = B\(dt*[Tx,Ty,Tz]*X(:,i));
end;

drNdu = sum(sensx,2);