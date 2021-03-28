function [ P, Tx, Ty, Tz] = dTrilinears(T,x,y,z,hx,hy,hz,bc)
% Rena 02/20/2017: allows mass to enter/leave through boundaries so mass is
% not necessarily conserved.
% Use linear interpolation to evaluate (approximate) the values of T in the
% grid (x,y,z)
% That = T(x,y,z) to be moved, necessary nly for the derivative
% Output: S, moves T, T^n+1 = S'*T^n;
%         Tx,Ty,Tz - derivatives of (S'*T^n) w.r.t. x, y and z coordinate
%%

if nargin < 8
    %bc = 'closed';
    bc = 'open';
end

[m,n,s] = size(x);

% Convert x and y to the coordinate system 1:m, 1:n, 1:s
% x = 1 + x/h1 + 1/2; y = 1+ y/h2 + 1/2; z = 1+ z/h3 + 1/2;
x = x./hx; y = y./hy; z = z./hz;
x = x + 1/2 ; y = y + 1/2 ; z = z + 1/2;

xd = floor(x(:)); xp  = x(:)-xd;
yd = floor(y(:)); yp = y(:)-yd;
zd = floor(z(:)); zp = z(:)-zd;
%{
ind1 = find(1<=xd & xd<m+1 & 1<=yd & yd<n+1 & 1<=zd & zd<s+1);
ind2 = find(1<=xd & xd<m   & 1<=yd & yd<n+1 & 1<=zd & zd<s+1);
ind3 = find(1<=xd & xd<m+1 & 1<=yd & yd<n   & 1<=zd & zd<s+1);
ind4 = find(1<=xd & xd<m   & 1<=yd & yd<n   & 1<=zd & zd<s+1);
ind5 = find(1<=xd & xd<m+1 & 1<=yd & yd<n+1 & 1<=zd & zd<s);
ind6 = find(1<=xd & xd<m   & 1<=yd & yd<n+1 & 1<=zd & zd<s);
ind7 = find(1<=xd & xd<m+1 & 1<=yd & yd<n   & 1<=zd & zd<s);
ind8 = find(1<=xd & xd<m   & 1<=yd & yd<n   & 1<=zd & zd<s);
%}
%ind = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);

switch bc
    case 'closed'
        %CLOSED VERSION dtri1 (use this one)
        %
        ind1 = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);
        ind2 = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);
        ind3 = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);
        ind4 = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);
        ind5 = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);
        ind6 = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);
        ind7 = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);
        ind8 = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);
        %}
    case 'open'
        % OPEN VERSION dtri3 (use this one)
        %
        ind1 = find(1<=xd & xd<m+1 & 1<=yd & yd<n+1 & 1<=zd & zd<s+1);
        ind2 = find(0<=xd & xd<m   & 1<=yd & yd<n+1 & 1<=zd & zd<s+1);
        ind3 = find(1<=xd & xd<m+1 & 0<=yd & yd<n   & 1<=zd & zd<s+1);
        ind4 = find(0<=xd & xd<m   & 0<=yd & yd<n   & 1<=zd & zd<s+1);
        ind5 = find(1<=xd & xd<m+1 & 1<=yd & yd<n+1 & 0<=zd & zd<s);
        ind6 = find(0<=xd & xd<m   & 1<=yd & yd<n+1 & 0<=zd & zd<s);
        ind7 = find(1<=xd & xd<m+1 & 0<=yd & yd<n   & 0<=zd & zd<s);
        ind8 = find(0<=xd & xd<m   & 0<=yd & yd<n   & 0<=zd & zd<s);
        %}
end

jki     =  xd(ind1) + (yd(ind1)-1)*m + (zd(ind1)-1)*m*n;
j1ki    =  xd(ind2) + 1+ (yd(ind2)-1)*m + (zd(ind2)-1)*m*n;
jk1i    =  xd(ind3) + yd(ind3)*m + (zd(ind3)-1)*m*n;
j1k1i   =  xd(ind4) + 1 + yd(ind4)*m + (zd(ind4)-1)*m*n;

jki1    =  xd(ind5) + (yd(ind5)-1)*m + zd(ind5)*m*n;
j1ki1   =  xd(ind6) + 1+ (yd(ind6)-1)*m + zd(ind6)*m*n;
jk1i1   =  xd(ind7) + yd(ind7)*m + zd(ind7)*m*n;
j1k1i1  =  xd(ind8) + 1 + yd(ind8)*m + zd(ind8)*m*n;

ii = [ind1;ind2;ind3;ind4;ind5;ind6;ind7;ind8];
jj = [jki; j1ki;   jk1i; j1k1i;   jki1; j1ki1;   jk1i1; j1k1i1];

ss = [(1-xp(ind1)).*(1-yp(ind1)).*(1-zp(ind1));...
    (xp(ind2)).*(1-yp(ind2)).*(1-zp(ind2));...
    (1-xp(ind3)).*(yp(ind3)).*(1-zp(ind3));...
    (xp(ind4)).*(yp(ind4)).*(1-zp(ind4));...
    (1-xp(ind5)).*(1-yp(ind5)).* (zp(ind5));...
    (xp(ind6)).*(1-yp(ind6)).* (zp(ind6)); ...
    (1-xp(ind7)).*(yp(ind7)).* (zp(ind7));...
    (xp(ind8)).*(yp(ind8)).* (zp(ind8))];
%%
%{
S = sparse(ii,jj,ss);
St = S';

P = sparse(n*m*s,n*m*s);
P(1:size(St,1),1:size(St,2)) = St;
P = squeeze(P(1:n*m*s,1:n*m*s));
%
%S = sparse(ii,jj,ss);
%St = S';
St=sparse(jj,ii,ss);

P = sparse(n*m*s,n*m*s);
P = St(1:n*m*s,1:n*m*s);
%}
%S = sparse(ii,jj,ss,n*m*s,n*m*s);
P=sparse(jj,ii,ss,n*m*s,n*m*s);
if nargout > 1
    
    %% Derivatives
    ind = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);
    
    jki     =  xd(ind) + (yd(ind)-1)*m + (zd(ind)-1)*m*n;
    j1ki    =  xd(ind) + 1+ (yd(ind)-1)*m + (zd(ind)-1)*m*n;
    jk1i    =  xd(ind) + yd(ind)*m + (zd(ind)-1)*m*n;
    j1k1i   =  xd(ind) + 1 + yd(ind)*m + (zd(ind)-1)*m*n;
    
    jki1    =  xd(ind) + (yd(ind)-1)*m + zd(ind)*m*n;
    j1ki1   =  xd(ind) + 1+ (yd(ind)-1)*m + zd(ind)*m*n;
    jk1i1   =  xd(ind) + yd(ind)*m + zd(ind)*m*n;
    j1k1i1  =  xd(ind) + 1 + yd(ind)*m + zd(ind)*m*n;
    
    A1 = sparse(jki,ind,T(ind),n*m*s,n*m*s);
    A2 = sparse(j1ki,ind,T(ind),n*m*s,n*m*s);
    A3 = sparse(jk1i,ind,T(ind),n*m*s,n*m*s);
    A4 = sparse(j1k1i,ind,T(ind),n*m*s,n*m*s);
    A5 = sparse(jki1,ind,T(ind),n*m*s,n*m*s);
    A6 = sparse(j1ki1,ind,T(ind),n*m*s,n*m*s);
    A7 = sparse(jk1i1,ind,T(ind),n*m*s,n*m*s);
    A8 = sparse(j1k1i1,ind,T(ind),n*m*s,n*m*s);
    
    v1dxz = zeros(n*m*s,1); v2dxz = zeros(n*m*s,1);
    v3dxz = zeros(n*m*s,1); v4dxz = zeros(n*m*s,1);
    v5dxdz = zeros(n*m*s,1); v6dxdz = zeros(n*m*s,1);
    v7dxdz = zeros(n*m*s,1); v8dxdz = zeros(n*m*s,1);
    
    v1dyz = zeros(n*m*s,1); v2dyz = zeros(n*m*s,1);
    v3dyz = zeros(n*m*s,1); v4dyz = zeros(n*m*s,1);
    v5dydz = zeros(n*m*s,1); v6dydz = zeros(n*m*s,1);
    v7dydz = zeros(n*m*s,1); v8dydz = zeros(n*m*s,1);
    
    v1dzz = zeros(n*m*s,1); v2dzz = zeros(n*m*s,1);
    v3dzz = zeros(n*m*s,1); v4dzz = zeros(n*m*s,1);
    v5dzdz = zeros(n*m*s,1); v6dzdz = zeros(n*m*s,1);
    v7dzdz = zeros(n*m*s,1); v8dzdz = zeros(n*m*s,1);
    
    v1dxz(ind)  = (-1)*(1-yp(ind)).*(1-zp(ind));
    v2dxz(ind)  = (1-yp(ind)).*(1-zp(ind));
    v3dxz(ind)  = (-1)*(yp(ind)).*(1-zp(ind));
    v4dxz(ind)  = (yp(ind)).*(1-zp(ind));
    v5dxdz(ind) = (-1)*(1-yp(ind)).* (zp(ind));
    v6dxdz(ind) = (1-yp(ind)).* (zp(ind));
    v7dxdz(ind) = (-1)*(yp(ind)).* (zp(ind));
    v8dxdz(ind) = (yp(ind)).* (zp(ind));
    
    
    Txx = A1*sdiag(v1dxz) + A2*sdiag(v2dxz)+ A3*sdiag(v3dxz) + A4*sdiag(v4dxz)+...
        A5*sdiag(v5dxdz) + A6*sdiag(v6dxdz)+ A7*sdiag(v7dxdz) + A8* sdiag(v8dxdz);
    
    v1dyz(ind)  = (1-xp(ind)).*(-1).*(1-zp(ind));
    v2dyz(ind)  = (xp(ind)).*(-1).* (1-zp(ind));
    v3dyz(ind)  = (1-xp(ind)).*(1-zp(ind));
    v4dyz(ind)  =  xp(ind).*(1-zp(ind));
    v5dydz(ind) = (1-xp(ind)).*(-1).*(zp(ind));
    v6dydz(ind) = (xp(ind)).*(-1).* (zp(ind));
    v7dydz(ind) = (1-xp(ind)).* (zp(ind));
    v8dydz(ind) =  xp(ind).*zp(ind);
    
    Tyy = A1*sdiag(v1dyz) + A2*sdiag(v2dyz)+ A3*sdiag(v3dyz) + A4*sdiag(v4dyz)+...
        A5*sdiag(v5dydz) + A6*sdiag(v6dydz)+ A7*sdiag(v7dydz) + A8* sdiag(v8dydz);
    
    v1dzz(ind) = (1-xp(ind)).* (1-yp(ind)).*(-1);
    v2dzz(ind) = xp(ind).* (1-yp(ind)).*(-1);
    v3dzz(ind) = (1-xp(ind)).*(yp(ind)).*(-1);
    v4dzz(ind) = (xp(ind)).*(yp(ind)).*(-1);
    v5dzdz(ind) = (1-xp(ind)).*(1-yp(ind));
    v6dzdz(ind) = xp(ind).*(1-yp(ind));
    v7dzdz(ind) = (1-xp(ind)).*yp(ind);
    v8dzdz(ind) = xp(ind).*yp(ind);
    
    Tzz = A1*sdiag(v1dzz) + A2*sdiag(v2dzz)+ A3*sdiag(v3dzz) + A4*sdiag(v4dzz)+...
        A5*sdiag(v5dzdz) + A6*sdiag(v6dzdz)+ A7*sdiag(v7dzdz) + A8* sdiag(v8dzdz);
    
    Tx = Txx/hx;
    Ty = Tyy/hy;
    Tz = Tzz/hz;
end;

end


%{
function [ P, Tx, Ty, Tz] = dTrilinears(T,x,y,z,hx,hy,hz)
% original version from Klara:
% Use linear interpolation to evaluate (approximate) the values of T in the
% grid (x,y,z)
% That = T(x,y,z) to be moved, necessary nly for the derivative
% Output: S, moves T, T^n=1 = S'*T^n;
%         Tx,Ty,Tz - derivatives of (S'*T^n) w.r.t. x, y and z coordinate
%%
[m,n,s] = size(x);

% Convert x and y to the coordinate system 1:m, 1:n, 1:s
% x = 1 + x/h1 + 1/2; y = 1+ y/h2 + 1/2; z = 1+ z/h3 + 1/2;
x = x./hx; y = y./hy; z = z./hz;
x = x + 1/2 ; y = y + 1/2 ; z = z + 1/2;

xd = floor(x(:)); xp  = x(:)-xd;
yd = floor(y(:)); yp = y(:)-yd;
zd = floor(z(:)); zp = z(:)-zd;
%%
ind = find(1<=xd & xd<m+1 & 1<=yd & yd<n+1 & 1<=zd & zd<s+1);
%ind = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);

jki     =  xd(ind) + (yd(ind)-1)*m + (zd(ind)-1)*m*n;
j1ki    =  xd(ind) + 1+ (yd(ind)-1)*m + (zd(ind)-1)*m*n;
jk1i    =  xd(ind) + yd(ind)*m + (zd(ind)-1)*m*n;
j1k1i   =  xd(ind) + 1 + yd(ind)*m + (zd(ind)-1)*m*n;

jki1    =  xd(ind) + (yd(ind)-1)*m + zd(ind)*m*n;
j1ki1   =  xd(ind) + 1+ (yd(ind)-1)*m + zd(ind)*m*n;
jk1i1   =  xd(ind) + yd(ind)*m + zd(ind)*m*n;
j1k1i1  =  xd(ind) + 1 + yd(ind)*m + zd(ind)*m*n;

ii = [ind;ind;ind;ind;ind;ind;ind;ind];
jj = [jki; j1ki;   jk1i; j1k1i;   jki1; j1ki1;   jk1i1; j1k1i1];

ss = [(1-xp(ind)).*(1-yp(ind)).*(1-zp(ind));...
    (xp(ind)).*(1-yp(ind)).*(1-zp(ind));...
    (1-xp(ind)).*(yp(ind)).*(1-zp(ind));...
    (xp(ind)).*(yp(ind)).*(1-zp(ind));...
    (1-xp(ind)).*(1-yp(ind)).* (zp(ind));...
    (xp(ind)).*(1-yp(ind)).* (zp(ind)); ...
    (1-xp(ind)).*(yp(ind)).* (zp(ind));...
    (xp(ind)).*(yp(ind)).* (zp(ind))];
%%
%{
S = sparse(ii,jj,ss);
St = S';

P = sparse(n*m*s,n*m*s);
P(1:size(St,1),1:size(St,2)) = St;
P = squeeze(P(1:n*m*s,1:n*m*s));
%}
%S = sparse(ii,jj,ss);
%St = S';
%St=sparse(jj,ii,ss);

%P = sparse(n*m*s,n*m*s);
%P = St(1:n*m*s,1:n*m*s);
%S = sparse(ii,jj,ss,n*m*s,n*m*s);
P=sparse(jj,ii,ss,n*m*s,n*m*s);
if nargout > 1
 
%% Derivatives
ind = find(1<=xd & xd<m & 1<=yd & yd<n & 1<=zd & zd<s);

jki     =  xd(ind) + (yd(ind)-1)*m + (zd(ind)-1)*m*n;
j1ki    =  xd(ind) + 1+ (yd(ind)-1)*m + (zd(ind)-1)*m*n;
jk1i    =  xd(ind) + yd(ind)*m + (zd(ind)-1)*m*n;
j1k1i   =  xd(ind) + 1 + yd(ind)*m + (zd(ind)-1)*m*n;

jki1    =  xd(ind) + (yd(ind)-1)*m + zd(ind)*m*n;
j1ki1   =  xd(ind) + 1+ (yd(ind)-1)*m + zd(ind)*m*n;
jk1i1   =  xd(ind) + yd(ind)*m + zd(ind)*m*n;
j1k1i1  =  xd(ind) + 1 + yd(ind)*m + zd(ind)*m*n;

A1 = sparse(jki,ind,T(ind),n*m*s,n*m*s);
A2 = sparse(j1ki,ind,T(ind),n*m*s,n*m*s);
A3 = sparse(jk1i,ind,T(ind),n*m*s,n*m*s);
A4 = sparse(j1k1i,ind,T(ind),n*m*s,n*m*s);
A5 = sparse(jki1,ind,T(ind),n*m*s,n*m*s);
A6 = sparse(j1ki1,ind,T(ind),n*m*s,n*m*s);
A7 = sparse(jk1i1,ind,T(ind),n*m*s,n*m*s);
A8 = sparse(j1k1i1,ind,T(ind),n*m*s,n*m*s);

v1dxz = zeros(n*m*s,1); v2dxz = zeros(n*m*s,1);
v3dxz = zeros(n*m*s,1); v4dxz = zeros(n*m*s,1);
v5dxdz = zeros(n*m*s,1); v6dxdz = zeros(n*m*s,1);
v7dxdz = zeros(n*m*s,1); v8dxdz = zeros(n*m*s,1);

v1dyz = zeros(n*m*s,1); v2dyz = zeros(n*m*s,1);
v3dyz = zeros(n*m*s,1); v4dyz = zeros(n*m*s,1);
v5dydz = zeros(n*m*s,1); v6dydz = zeros(n*m*s,1);
v7dydz = zeros(n*m*s,1); v8dydz = zeros(n*m*s,1);

v1dzz = zeros(n*m*s,1); v2dzz = zeros(n*m*s,1);
v3dzz = zeros(n*m*s,1); v4dzz = zeros(n*m*s,1);
v5dzdz = zeros(n*m*s,1); v6dzdz = zeros(n*m*s,1);
v7dzdz = zeros(n*m*s,1); v8dzdz = zeros(n*m*s,1);

v1dxz(ind)  = (-1)*(1-yp(ind)).*(1-zp(ind));
v2dxz(ind)  = (1-yp(ind)).*(1-zp(ind));
v3dxz(ind)  = (-1)*(yp(ind)).*(1-zp(ind));
v4dxz(ind)  = (yp(ind)).*(1-zp(ind));
v5dxdz(ind) = (-1)*(1-yp(ind)).* (zp(ind));
v6dxdz(ind) = (1-yp(ind)).* (zp(ind));
v7dxdz(ind) = (-1)*(yp(ind)).* (zp(ind));
v8dxdz(ind) = (yp(ind)).* (zp(ind));


Txx = A1*sdiag(v1dxz) + A2*sdiag(v2dxz)+ A3*sdiag(v3dxz) + A4*sdiag(v4dxz)+...
    A5*sdiag(v5dxdz) + A6*sdiag(v6dxdz)+ A7*sdiag(v7dxdz) + A8* sdiag(v8dxdz);

v1dyz(ind)  = (1-xp(ind)).*(-1).*(1-zp(ind));
v2dyz(ind)  = (xp(ind)).*(-1).* (1-zp(ind));
v3dyz(ind)  = (1-xp(ind)).*(1-zp(ind));
v4dyz(ind)  =  xp(ind).*(1-zp(ind));
v5dydz(ind) = (1-xp(ind)).*(-1).*(zp(ind));
v6dydz(ind) = (xp(ind)).*(-1).* (zp(ind));
v7dydz(ind) = (1-xp(ind)).* (zp(ind));
v8dydz(ind) =  xp(ind).*zp(ind);

Tyy = A1*sdiag(v1dyz) + A2*sdiag(v2dyz)+ A3*sdiag(v3dyz) + A4*sdiag(v4dyz)+...
    A5*sdiag(v5dydz) + A6*sdiag(v6dydz)+ A7*sdiag(v7dydz) + A8* sdiag(v8dydz);

v1dzz(ind) = (1-xp(ind)).* (1-yp(ind)).*(-1);
v2dzz(ind) = xp(ind).* (1-yp(ind)).*(-1);
v3dzz(ind) = (1-xp(ind)).*(yp(ind)).*(-1);
v4dzz(ind) = (xp(ind)).*(yp(ind)).*(-1);
v5dzdz(ind) = (1-xp(ind)).*(1-yp(ind));
v6dzdz(ind) = xp(ind).*(1-yp(ind));
v7dzdz(ind) = (1-xp(ind)).*yp(ind);
v8dzdz(ind) = xp(ind).*yp(ind);

Tzz = A1*sdiag(v1dzz) + A2*sdiag(v2dzz)+ A3*sdiag(v3dzz) + A4*sdiag(v4dzz)+...
    A5*sdiag(v5dzdz) + A6*sdiag(v6dzdz)+ A7*sdiag(v7dzdz) + A8* sdiag(v8dzdz);

Tx = Txx/hx;
Ty = Tyy/hy;
Tz = Tzz/hz;
end;

end






%}

