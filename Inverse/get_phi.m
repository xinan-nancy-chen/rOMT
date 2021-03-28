%function[phi,mk,phi0,phiN] =  get_phi(rho0,u,nt,dt,par)
function[phi,mk,phiN,rho,Ru] =  get_phi(rho0,u,nt,dt,par)


rho   = advecDiff(rho0,u,nt,dt,par);

mk    = MKdist(u,rho,nt,dt, par); %=hd*dt*rho'*||v||^2

phiN  = 0.5*norm(rho(:,end) - par.drhoN)^2;
%%%%MASKED??????%%%%%
%% smoothing deformation field
%
Ru = 0;

if par.gamma~=0
    uvec=vec2mat(u(:),3*nt);
    %GTG=par.Grad'*par.Grad;
    for ii = 1:3*nt
        %Ru=Ru+0.5*par.hd*dt*(uvec(:,ii)'*GTG*uvec(:,ii));
        Ru=Ru+0.5*par.hd*dt*(norm(par.Grad*uvec(:,ii))^2); %Rgradu+0.5*(uvec(:,ii)'*uvec(:,ii));
    end
    %}
end

%phi   = par.beta*mk + par.alpha*phi0 + par.alpha*phiN + par.gamma.*Ru;
phi   = par.beta*mk + phiN + par.gamma.*Ru;
%phi   = par.beta*mk  + 0.5*norm(rho0 - par.drho0)...
%        + 0.5*norm(rho(:,end) - par.drhoN);

