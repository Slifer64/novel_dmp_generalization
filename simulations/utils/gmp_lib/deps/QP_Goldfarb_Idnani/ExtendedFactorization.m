function [C,Error]=ExtendedFactorization(C,tol)
% [C,Error]=ExtendedFactorization(C,tol)
% Cholesky like decomposition of possibly semi definite matrix
%
%   The code factorizes a symmetric C by an triangular decompo-
%   sition, C = U'*U, where U is an upper triagonal matrix. In case of a
%   semi-definite matrix, or almost negative definite, the orginal C matrix
%   is modified such that its smaller eigenvalue is larger or equal to tol
%   If the given C has its smaller eigenvalue less than tol, C is considered
%   indefinite, no U is computed and Error is set to true.
%
%  Input parameters:
%   C(n,n)    :  positive definite matrix which is to be decomposed. note that
%                     only the upper triangular part of C including diagonal is read.
%   tol          : positive definite tolerance,(default ~1.0e-12).
%
% Output parameters
%   C(n,n)    :  Contains the upper triangular factor U of C in the
%                     upper triangular part including diagonal.
%   Error      :
%                   false - the decomposition is successfully computed.
%                   true  - C contains the original matrix
%
%  abc.consultant@wanadoo.fr copyright 2018
%  Author : Alain Barraud
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2, tol=eps^.75;end
n=size(C,1);
[V,D]=eig(C);% C=V*D*V';
d=diag(D);
if any(d<-tol)
    Error=true;
else
    Error=false;
    d=sqrt(max(d,0));
    V=V';%==> C=V'*D*V;
    for i=1:n, V(i,:)=d(i)*V(i,:);end
    C=triu(qr(V));
    L=diag(C)<0;C(L,:)=-C(L,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%