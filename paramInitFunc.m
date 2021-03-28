function par_out = paramInitFunc(n,nt,dt,sigma,add_source,gamma,beta,niter_pcg)
% this function returns the par structure

%  defines the models size, parameters
%  BC values if needed

%% Domain size, uniform grid
%n = dpar.true_size';
n1 = n(1);
n2 = n(2); 
n3 = n(3);

h1 = ones(n1,1);
h2 = ones(n2,1);
h3 = ones(n3,1);

h  = [h1;h2;h3];

[Xc, Yc, Zc] = getCellCenteredGrid(h1,h2, h3);

%% Matrix operators
BC                 = {'ccn' 'ccn' 'ccn'};
Grad               = getCellCenteredGradMatrix(BC,h1,h2, h3);                

%% add source:
if add_source
    qex = zeros(n');
    qex(qexi) = si;
    par.qex = sparse(qex(:));
end

%% Put it in a structure
par.add_source = add_source;
par.n     = n;
par.h1    = h1; par.h2 = h2; par.h3 = h3;
par.hd    = h1(1)*h2(1)*h3(1); %% NOTE: assuming all entries in h1 are the same, and the same for h2,h3
par.dt    = dt;
par.nt    = nt;
par.sigma = sigma;
par.Grad  = Grad; 
par.Xc    = Xc; par.Yc = Yc; par.Zc = Zc;

%% Inversion parameters
par.gamma     = gamma;% regularization term to smoothe the deformation field
par.alpha     = 0.1;%1;%0.1;
par.beta      = beta;
par.minUpdate = 1;
par.maxCGiter = 10;
par.niter_pcg = niter_pcg;
par.maxUiter  = 6;

par_out = par;
end
