function[mk] = MKdist(u,rho,nt,dt, par)
%[I] = MKdist(u,rho,n,nt,h,dt)
%
% build averaging matrix
n   = par.n;

A   = kron(ones(1,3),speye(prod(n)));

U   = reshape(u,3*prod(n),nt);
R   = reshape(rho,prod(n),nt);

mk = 0;


for i=1:nt
    mk = mk + par.hd*dt*R(:,i)'*A*(U(:,i).*U(:,i));
end

%mk = dt*R(:,end)'*A*(U(:,end).*U(:,end));
%}
%% OR EQUIVALENTLY: (DOUBLE CHECK EQUIV WHEN nt > 1)
%{
%function[mk] = MKdist(u,rho,nt,dt, par)
n   = par.n;
A   = kron(speye(nt),kron(ones(1,3),speye(prod(n))));

mk  = par.hd*dt*rho'*A*(u.*u);
%}