function [mk,dmkdu] = get_dmkdu(rho,rho_0,u,nt,dt,par)
%% Sensitivities of mk w.r.t u
%  x - vector size (nt*3*prod(n),1)

n   = par.n;
A   = kron(ones(1,3),speye(prod(n)));

U   = reshape(u,3*prod(n),nt);
R   = reshape(rho,prod(n),nt);

mk   = 0;
dmku = zeros(3*prod(n),nt);


for i = 1:nt
    mk = mk + dt*R(:,i)'*A*(U(:,i).*U(:,i));
   
    dmku(:,i) = dmku(:,i) + 2*dt*(R(:,i)'*A*diag(sparse(U(:,i))))';
    
    dmku(:,1:i) = dmku(:,1:i) + reshape(dt*get_drNduT(rho_0,U(:,1:i),i,dt,par,A*(U(:,i).*U(:,i))),3*prod(n),i);
                   
end

dmkdu = dmku(:);