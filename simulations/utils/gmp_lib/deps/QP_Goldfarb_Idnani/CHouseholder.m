function [v,u,b]=CHouseholder(v,Lz,piv)
%  [v,u,b]=CHouseholder(v,Lz,piv)
% compute an elementary Householder transformation Q
% such that the Lz list of components of the new v are  zeros,  piv is the pivot index
% component of v. The other components of v are unchanged.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT
% v         : v the input vector
% Lz       : the indices list of v components which must be zeroed
% piv      : the v pivot index
%
% OUTPUT
% v         : the transformed vector, with v(Lz)=0
% u         : such that Q = I - 2*u*u'/(u'*u)
% b         : a pre computed value used to speed up Q*w or w*Q products
%              where w is a column or row vector
%
% The Linpack algorithm is used, v is real or complex
%
% Author : Alain Barraud, copyright abc.consultant@wanadoo.fr 2003-2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(size(v));%initialization and memory allocation
vp=v(piv);rv1=real(vp);%the real pivot part
if isrow(v), s=norm([v(piv),v(Lz)]);else, s=norm([v(piv);v(Lz)]);end%s computation
if rv1~= 0, s = s*sign(rv1);end
u(Lz) = v(Lz);%zeroed part of v
up = vp + s;b=s*up;%compute u(piv)
u(piv)=up;%now u is updated
if b~=0%Non trivial case Q~=I
    v(piv) = -s;v(Lz)=0;b = -1/b;%update new v
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%