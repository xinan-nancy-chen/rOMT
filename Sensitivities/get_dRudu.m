function dRdu=get_dRudu(u,nt,par)
% -----------------------------
% Created by:       Rena Elkin
% Last modified on: 04/13/2017
% -----------------------------
%This function returns the derivative of the regularization term
%R(u)=||grad(u)||^2

%u=u(:);
G=-par.Grad'*par.Grad;
U=vec2mat(u(:),3*nt);
dRdu=G*U;
dRdu=dRdu(:)';

%% dRdu = G*u(:);
