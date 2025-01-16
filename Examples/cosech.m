function mat = cosech(M,Val,Minv)
% mat = cosech(M,Val,Minv)
% COSECH: Cosecante hiperbolica de una matriz X.
% M     : Matriz de eigenvectores de X.
% Minv  : Matriz inversa de M.
% Val   : Vector de valores caracteristicos de X.
%
arg = exp(-Val);
dd = 2.0*arg./(1 - arg.*arg);
mat = M*diag(dd)*Minv;