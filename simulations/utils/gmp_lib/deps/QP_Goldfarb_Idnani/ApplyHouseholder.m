function M=ApplyHouseholder(M,u,b,side)
% M=ApplyHouseholder(M,u,b,side)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Apply an elementary Householder transformation Q from the RIGHT or from
% the LEFT, according to side,  to a set of vectors packed as a matrix M
% u and b are previously computed by CHouseholder
%
% INPUT
% M        : the input matrix, size mxn, considered as n columns ('Left') or m rows ('Right')
% u         : such that Q = I - 2*u*u'/(u'*u),  mx1 ('Left')  or size nx1 ('Right')
% b        : a pre computed value used to speed up Q*w(:)  product
% side   : specifies the side upon which Q is applied, recognised values :
%            'L' stands for 'Left' multiplication, default value
%            'R' stands for 'Right' multiplication
%
% OUTPUT
% M  : the transformed matrix  Q*M ('Left') or M*Q ('Right')
%
% Author : Alain Barraud, copyright abc.consultant@wanadoo.fr 2003-2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%force upper case
if nargin==3, side='L';else, side=upper(side(1));end
I = find(u~=0);
if b==0||isempty(I); return;end
v=u(I);%only non trivial part of u is used
switch side
    case 'R'
        tau=v*b;
        tau=M(:,I)*tau;%intermediate row vector
        M(:,I) =M(:,I) + tau*v';%transformed M
    case 'L'
        tau=v'*b;tau=tau*M(I,:);%intermediate row vector
        M(I,:) =M(I,:) + v*tau;%transformed M
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%