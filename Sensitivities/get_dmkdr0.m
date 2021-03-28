function[mk,dmk, d2mk] = get_dmkdr0(u,rho,rho0,nt,dt, par)

n   = par.n;
A   = kron(ones(1,3),speye(prod(n)));

U   = reshape(u,3*prod(n),nt);
R   = reshape(rho,prod(n),nt);

mk  = 0;
dmk = 0;

for i=1:nt
    %mk  = mk + dt*R(:,i)'*A*(U(:,i).*U(:,i));
    dmk = dmk + par.hd*dt*get_drNdr0T(rho0,U(:,1:i),i,dt,par,A*(U(:,i).*U(:,i)));   
end

d2mk = [];


%% or alternatively:
%{
function dmk = get_dmkdr0(u,nt,dt,par)
n                  = par.n;
A                  = kron(speye(nt),kron(ones(1,3),speye(prod(n))));
M                  = speye(prod(n));
Mdis               = -par.sigma*par.Grad'*par.Grad;
B                  = inv(M - dt*Mdis);
dmk                = zeros(prod(n),1);
x                  = vec2mat(A*(u(:).*u(:)),nt);
u                  = reshape(u(:),3*prod(n),[]);
for i = 1:nt
    U1 = reshape(u(1:prod(n),i),n');
    U2 = reshape(u(prod(n)+1:2*prod(n),i),n');
    U3 = reshape(u(2*prod(n)+1:end,i),n');
    S  = dTrilinears(ones(prod(n),1),par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3, ...
                     par.h1(1),par.h2(1),par.h3(1));
    M   = S'*B*M;
    dmk = dmk + M*x(:,i);
end
dmk = par.hd*dt*dmk;
%}