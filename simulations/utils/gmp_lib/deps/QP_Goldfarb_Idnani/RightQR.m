function varargout=RightQR(varargin)
%
% [R,Q]=RightQR(A)
% Computes an orthogonal factorization Q such that A*Q = R is lower triangular
% by means of householder transformations
%
% when a third output argument is added a row permutation is done such that
% diagonal elements of R are sorted in absolute decreasing values
% [R,Q,e]=RightQR(A),   A(:,e)*Q = R
%
%  Remark .- Real and complex matrix ar supported
%
% Author : Alain Barraud, copyright abc.consultant@wanadoo.fr 2003-2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
narg=nargin;
if narg==0,error('RightQR requires at  least one argument');end
A=varargin{1};
if nargout==3, Pivoting=true;else, Pivoting=false;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialisation and memory allocation
[m,n] = size(A);%A dimension
if Pivoting, e=(1:m);end
r=min(m,n);
Q=eye(n);%memory allocation
%begin factorization
for k = 1:r%loop on the A rows
    if Pivoting
        J=(k:n);%the set of used columns indices for pivoting
        I=(k:m);%the set of used rows indices for pivoting
        tmp=sum(A(I,J).*conj(A(I,J)),2);[tmpx,piv]=max(tmp);
        piv=k-1+piv;%piv is the next row to be used
        if piv~=k, tmp=e(k);e(k)=e(piv);e(piv)=tmp;tmp=A(k,:);A(k,:)=A(piv,:);A(piv,:)=tmp;end
    end
    [A(k,:),uk,betak]=CHouseholder(A(k,:),k+1:n,k);%compute the kth Householder transformation
    if Pivoting
        if k==1
            Ref=abs(A(1,1));
        elseif tmpx<n*eps*Ref
            A(k:m,:)=0;
            break
        end
    end
    %apply this transformation to rows A(I,:) if any
    if k<m, I=(k+1:m);A(I,:)=ApplyHouseholder(A(I,:),uk(:),betak,'R');end
    Q=ApplyHouseholder(Q,uk(:),betak,'R');
end
varargout{1}=A(1:r,:);varargout{2}=Q;
if nargout==3, varargout{3}=e;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%