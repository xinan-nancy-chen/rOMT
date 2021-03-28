function [u,phi,g] = GNblock_u(Rho_i,u,nt,dt,par,tag_str)
%%

if nargin < 6
    tag_str = '';
end

phi          = get_phi(Rho_i,u,nt,dt,par);
A            = kron(ones(1,3),speye(prod(par.n)));
Abig         = kron(speye(nt),A);
flag         = 0;
U            = reshape(u,3*prod(par.n),[]);
dmk2         = zeros(3*prod(par.n),nt);

%% Loop
for i = 1:par.maxUiter
    Rho      = advecDiff(Rho_i,u,nt,dt,par);
    for      j = 1:nt
             dmk2(:,1:j) = dmk2(:,1:j) + reshape(par.hd*dt*get_drNduT(Rho_i,U(:,1:j),j,dt,par,A*(U(:,j).*U(:,j))),3*prod(par.n),j);
    end
    
    g        = (par.beta*2*par.hd*dt*Rho(:)'*Abig*sdiag(u(:)))' + par.beta*dmk2(:) + ...
                get_drNduT(Rho_i,u,nt,dt,par,Rho(:,end) - par.drhoN) + par.gamma*dt*par.hd*get_dRudu(u,nt,par)';
    
    fprintf('%3d.%d\t      %3.2e \t     ||g|| = %3.2e       %s\n',i,0,phi,norm(g),tag_str);
           
    H        = @(x) par.beta*2*dt*par.hd*(Rho(:)'*Abig*sdiag(x))' + ...
               get_drNduT(Rho_i,u,nt,dt,par,get_drNdu(Rho_i,u,nt,dt,par,x))+...
               par.gamma*dt*par.hd.*kron(speye(nt*3),par.Grad'*par.Grad)*x;      
 
    
    [s,pcgflag,relres,iter]    = pcg(H,-g,0.01,par.niter_pcg);
    
    if pcgflag ~= 0
      warning('MATLAB:pcgExitFlag','Warning: GNblock_u.m >>> iter %d, while finding s, pcg exit flag = %d \nrelres = %3.2e, iter = %d, %s',i,pcgflag,relres,iter,tag_str)
    end
    
    muls     = 0.7; lsiter = 1;
    while 1
        ut   = u(:) + muls*s;
        
        phit = get_phi(Rho_i,ut,nt,dt,par);
        
        fprintf('%3d.%d\t      %3.2e \t     phit  = %3.2e        %s\n',i,lsiter,phi,phit,tag_str);
        
        % test for line search termination
        if phit < phi + 1e-8*s'*g
            break;                      %breaks while loop entirely (and goes to next statement after end of while loop)
        end
        muls = muls/2; lsiter = lsiter+1;
        
        % fail if lsiter is too large
        if lsiter > 4
            fprintf('LSB\n');
            ut = u;     
            flag = 1;    
            return;                     % returns and exits function
        end
    end
    
    if flag
        return; 
    end                % returns and exits function
    u   = ut;
    phi = phit;
    
end

end
