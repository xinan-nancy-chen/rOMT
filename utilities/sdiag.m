function [S] = sdiag(d)
% Returns d on the main diag of a sparse matrix.
S = spdiags(d(:),0,numel(d),numel(d));