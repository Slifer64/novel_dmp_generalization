function [x,fval,exitflag,output,lambda]=qpGoldfarbIdnani(H,f,Ai,bi,Ae,be,lb,ub,option)
% [x,fmin,errflag,output,lambda]=QP(H,f,Ai,bi,Ae,be,bi,lb,ub,option)
% [x,fmin,errflag,output,lambda]=QP(H,f,Ai,bi,Ae,be,lb,ub)
% [x,fmin,errflag,output,lambda]=QP(H,f,Ai,bi,Ae,be,lb)
% [x,fmin,errflag,output,lambda]=QP(H,f,Ai,bi,Ae,be)
% [x,fmin,errflag,output,lambda]=QP(H,f,Ai,bi)
% [x,fmin,errflag,output,lambda]=QP(H,f)
%
% Solves dense General Convex Quadratic Program :
% min 0.5*x'*H*x + f'*x ,  such that Ae*x=be; Ai*x<=bi;  lb<=x<=ub
%  x
%  with the following hypothesis :H'=H>=0 and rank([H;Ae])=n
%  i.e. reduced Hessian > 0 when Ae is not empty
%                       Hessian H > 0 when Ae is empty
%
% Input
% H          : cost H term, H'=H>=0 matrix, size nxn
% f           : cost f term , size nx1
% Ai        : The inequality constraint matrix, size mxn, m>=0
% bi         : The right hand side, size mx1
% Ae        : The equality constraint matrix, size pxn, p>=0
% be         : The right hand side, size px1
% lb         : defines the lower bounds, lb<=x set lb(j)=-inf for any lower unbounded x(j); size nx1
%                when unspecified, the problem is considered as lower unbounded
% ub         : defines the upper bounds, x<=ub; set ub(j)=+inf for any upper unbounded x(j); size nx1
%                when unspecified, the problem is considered as unbounded
%                In any case ub>=lb, fixed variables are supported.
% option  : a structure with 2 fields
%                 .maxiter controls the maximum number of iterations, (default 10*(n+m+p);
%                 .tol tolerance for infeasibility and linear independance,  (default 1e-12)
%
% Ouput
%  x            : the  solution, size nx1
% fval         : the minimum found
% exitflag  : exit status
%                 1 - Proven Optimal
%                 0 - Maximum Iterations reached
%                -1 - Infeasible Linear Equality
%                -2 - Infeasible Linear Inequality
%                -3 - Reduced Hessain singular
%                -4 - Non ConvexQP Problem
% output    : a structure with 4 fields
%                 .status strings indicating what happens
%                 .iter the number of iterations to achieve convergence
%                 .drop the number of removed constraints
%                  .tcpu the cpu time
% lambda  : the Lagrange multipliers a structure with fields .lower, .upper, .eqlin, .ineqlin
%                 lambda>0 for an active inequality or bound constraint
%                 lambda=0 for inactive inequality or bound constraint
%                 lambda~=0 for equality constraints
%
%   Author Alain Barraud
%   abc.consultant Copyright 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2,n=5;f=-ones(n,1);end;if nargin<1,H=eye(n);end%to avoid error
n=length(f);f=f(:);%problem size
%set default arguments
if nargin<8, ub=inf(n,1);end
if nargin<7, lb=-inf(n,1);end
if nargin<6, Ae=[];be=[];end
if nargin<4, Ai=[];bi=[];end
%get constraints sizes
p=length(be);%number of equality constraints
m=length(bi);%number of inequality constraints
%check parameters
if length(ub)~=n||length(lb)~=n
    error(['lb and ub length must be ',num2str(n)]);
end
if p>0
    [pp,nn]=size(Ae);
    if pp~=p||nn~=n
        error(['Matrix Ae must have ',num2str(p),' rows and ',numstr(n),' columns']);
    else
        be=be(:);
    end
end
if m>0
    [mm,nn]=size(Ai);
    if mm~=m||nn~=n
        error(['Matrix Ai must have ',num2str(m),' rows and ',numstr(n),' columns']);
    else
        bi=bi(:);
    end
end
ub=ub(:);lb=lb(:);ub=max(ub,lb);
%set default options values
Maxiter=10*(n+m+p);tol=1e-12;
if nargin<9
    option.maxiter=Maxiter;option.tol=tol;
else
    if ~isfield(option,'maxiter'), option.maxiter=Maxiter;end
    if ~isfield(option,'tol'), option.tol=tol;end
end
tic;
%%%
if p>0%Reduce the original QP problem 
    [R,fr,Air,bir,y1,Q,Es,status,Flag]=ReduceQP(H,f,Ae,be,Ai,bi,lb,ub,tol);
    exitflag=Flag;
else
    [Es,Error]=ExtendedFactorization(H,option.tol);
    Flag=nan;
    if Error
        output.status='H is not positive definite';exitflag=-3;
        output.iter=[];output.drop=[];
    else
        exitflag=1;
        R=Es;fr=f;Air=Ai;bir=bi;y1=[];Q=eye(n);
    end
end
if exitflag<0
    x=[];fval=[];lambda=[];output.status=status;
else%solves a inequality contraints strictly convex QP
    mr=length(bir);%number of inequality constraints
    if mr>0%call Goldfard and Idnani algorithm
        [y2,~,exitflag,output,u2]=GoldfarfIdnani(R,fr,Air,bir,option);
    elseif p==0%unconstrained explicit solution
        y2=-H\f;u2=[];output.status='OK';output.iter=[];output.drop=[];
    else
        y2=[];u2=[];%only equality constrained problem
    end
    x=Q*[y1;y2];%full solution
    if Flag==2;output.status{2}=status;end
    fval=0.5*sum((Es*x).^2)+f'*x;%optimum value
    Lu=isfinite(ub);Ll=isfinite(lb);
    nu=sum(Lu);
    lambda.ineqlin=u2(1:m);
    lambda.upper=zeros(n,1);lambda.upper(Lu)=u2(m+1:m+nu);
    lambda.lower=zeros(n,1);lambda.lower(Ll)=u2(m+nu+1:mr);
    if p>0
        lambda.eqlin=-Ae'\(H*x+f+Ai'*lambda.ineqlin+lambda.upper-lambda.lower);
    else
        lambda.eqlin=[];
    end
    %check for fixed variables
    fixed=lb==ub;
    gx=H(fixed,:)*x+f(fixed);
    Lx=find(fixed);nx=length(Lx);
    %among the two equal bounds only one is active!!
    for k=1:nx
        val=max(lambda.lower(Lx(k)),lambda.upper(Lx(k)));
        %lower bound is active when gradient > 0
        if gx(k)>0 , lambda.lower(Lx(k))=val;lambda.upper(Lx(k))=0;end
        %upper bound is active when gradient < 0
        if gx(k)<0 , lambda.upper(Lx(k))=val;lambda.lower(Lx(k))=0;end
    end
end
output.tcpu=toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%