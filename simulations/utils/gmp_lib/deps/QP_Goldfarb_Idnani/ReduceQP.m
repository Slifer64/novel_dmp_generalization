function [R,fr,Air,bir,y1,Q,Es,Status,Flag]=ReduceQP(H,f,Ae,be,Ai,bi,lb,ub,tol)
% [Hr,fr,Air,bir,fvr,y1,Q,,Es,Status,Flag]=ReduceQP(H,f,Ae,be,Ai,bi,tol)
% Given the following general Convex QP problem :
% min 0.5*x'*H*x + f'*x ,  such that Ae*x=be; Ai*x<=bi;  lb<=x<=ub
%   x
% with H'=H>=0 and [Ae;H] full column rank
% ReduceQP computes a reduced strictly convex QP problem
% without any equality constraints :
% min 0.5*y2'*Hr*y2 + fr'*y2 ,  such that Air*y2<=bir;
%   y2
% and y1 such that the solution of the original problem is x=Q*[y1;y2]
% where y2 is the solution of the reduced problem of size (n-p).
%
% Note [p,n]=size(Ae), Q=[Q1,Q2] such that Hr=Q2'*H*Q2>0, Ae*Q2=0;
% Q1 is nxp and Q2 is nx(n-p)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(f);Flag=1;Status='OK';
if nargin<9, tol=16*n*eps;end
if nargin<8, ub=inf(n,1);end
if nargin<7, lb=-inf(n,1);end
if nargin<6, Ai=[];bi=[];end
if nargin<4, error('ReduceQP requires at least 4 parameters');end
[E,Error]=ExtendedFactorization(H,tol);Es=E;
if ~Error
    p=length(be);
    %Null space decomposition
    [R,Q,e]=RightQR(Ae);
    ref=abs(R(1,1))*tol;r=1;
    while r<p&&abs(R(r+1,r+1))>ref;r=r+1;end
    L1=1:r;L2=r+1:n;
    y1=R(L1,L1)\be(e(L1));%particular solution of equality constraints
    if r<p%check rank deficiency
        err=be-Ae*Q(:,L1)*y1;
        if norm(err)<tol*norm(be)%however feasible
            Status='Rank deficient Feasible Equality Constraints';
            Flag=2;
        else
            Flag=-1;
            Status='Infeasible Equality constraints';
           Air=[];bir=[];R=[];fr=[];y1=[];Q=[];Es=[];return;
        end
    end
    E2=E*Q(:,L2);%reduced hessian factorization
    [R,~,er]=LeftQR(E2);
    ref=abs(R(1,1))*tol;r=1;
    while r<n-p&&abs(R(r+1,r+1))>ref;r=r+1;end
    if r<n-p
        Flag=-3;
        Status='Reduce Hessian is singular';
        Air=[];bir=[];R=[];fr=[];y1=[];Q=[];Es=[];return;
    end
    Eye=eye(n);%Full Identity matrix
    R(:,er)=R;
    E1=E*Q(:,L1);
    fr=E2'*E1*y1+Q(:,L2)'*f;%reduced linear term
    %add bounds if any
    Lu=isfinite(ub);Ll=isfinite(lb);
    Ai=[Ai;Eye(Lu,:);-Eye(Ll,:)];bi=[bi;ub(Lu);-lb(Ll)];
    if ~isempty(Ai)
        Ai=Ai*Q;
        Air=Ai(:,L2);bir=bi-Ai(:,L1)*y1;
    else
        Air=[];bir=[];
    end
else 
    Flag=-4;
    Status='The QP problem is not convex';
    Air=[];bir=[];R=[];fr=[];y1=[];Q=[];Es=[];return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%