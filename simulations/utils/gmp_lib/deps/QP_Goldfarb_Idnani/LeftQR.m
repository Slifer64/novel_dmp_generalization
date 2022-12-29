function varargout=LeftQR(varargin)
%
% [R,Q]=LeftQR(A)
% Computes an orthogonal factorization Q such that Q*A = R is upper triangular
% by means of householder transformations
%
% when a third output argument is added a column permutation is done such that
% diagonal elements of R are sorted in absolute decreasing values
% [R,Q,e]=LeftQR(A),   Q*A(:,e)=R
%
% when a second input arguments is present, Q is not computed explicitly.
% orthogonal transformations are directly applied to B. 
% on output B contains the transformed B.
% [R,B]=LeftQR(A,B)
% [R,B,e]=LeftQR(A,B)
%
%  Remark .- Real and complex matrix ar supported
%
% Author : Alain Barraud, copyright abc.consultant@wanadoo.fr 2003-2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
narg=nargin;
if narg==0,error('LeftQR requires at  least one argument');end
A=varargin{1};
if narg==2
   B=varargin{2};
else
    B=[];
end
if nargout==3, Pivoting=true;else, Pivoting=false;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialisation and memory allocation
[m,n] = size(A);%A dimension
if Pivoting, e=(1:n);end
r=min(m,n);
if nargin==1
    Q=eye(m);%memory allocation
end
%begin factorization
for k = 1:r%loop on the A columns
    if Pivoting%do columns permutation (rank revealing)
        I=(k:m);%the set of used row indices for pivoting
        J=(k:n);%the set of used columns indices for pivoting
        tmp=sum(A(I,J).*conj(A(I,J)),1);[tmpx,piv]=max(tmp);
        piv=k-1+piv;%piv is the next column to be used
        if piv~=k, tmp=e(k);e(k)=e(piv);e(piv)=tmp;tmp=A(:,k);A(:,k)=A(:,piv);A(:,piv)=tmp;end
    end
    [A(:,k),uk,betak]=CHouseholder(A(:,k),(k+1:m),k);%compute the kth Householder transformation
    if Pivoting
        if k==1
            Ref=abs(A(1,1));
        elseif tmpx<n*eps*Ref
            A(k:m,:)=0;
            break
        end
    end
    %apply this transformation to the following A(:,J) columns if any
    if k<n,J=(k+1:n);A(:,J)=ApplyHouseholder(A(:,J),uk,betak,'L');end
    if nargin==1%Compute Q explicitly
        Q=ApplyHouseholder(Q,uk,betak,'L');
    elseif ~isempty(B)%apply Q to B
        B=ApplyHouseholder(B,uk,betak,'L');
    end
end
varargout{1}=A(1:r,:);
if nargin==2
    varargout{2}=B;
else
    varargout{2}=Q;
end
if nargout==3
    varargout{3}=e;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%