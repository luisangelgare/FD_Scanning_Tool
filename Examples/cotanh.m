function mat = cotanh(M,Val,Minv)
% mat = cotanh(M,Val,Minv)
% COTANH: Cotangente hiperbolica de una matriz.
% M     : Matriz de eigenvectores.
% Minv  : Matriz inversa de M.
% Val   : Vector de valores caracteristicos.
%
arg = exp(-2.0*Val);
dd = ( 1.0 + arg )./( 1.0 - arg );
mat = M*diag(dd)*Minv;