function [x,fval,exitflag,output,u]=GoldfarfIdnani(R,f,Ai,bi,option)
%
% Solves the following strictly Convex QP problem
% using Goldfarb Idnani algorithm
% min 0.5*x'*(R'*R)*x + f'*x ,  such that Ai*x<=bi;
% with H'=H=R'*R>0,  R is a square regular factorization of H (not necessarely Choleky)
%
% R,f,Ai,bi : Problem data 
% option    : see calling program
%
% x,fval      : problem solution and minimun
% exitflag   : 1 optimum found
%                   0 maximum iteration reached
%                  -2 Infeasible Inequality Constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=length(bi);n=length(f);
J=inv(R);
x=-R'\f;x=J*x;%#ok
fval=-0.5*norm(R*x)^2;
%initialisation
iq=1;R=[];exitflag=1;output.drop=0;
FullStep= true;W=[];
ListIn=(1:m);opts.UT=true;iter=0;drop=0;u=zeros(m,1);
while iter<=option.maxiter%main finite iterative loop
    nact=setxor(ListIn,W);
    %%% STEP1: check for violated constraints  
    if FullStep
        ctr = SelectConstraint(x,Ai,bi,nact,option.tol);
        %all constraints are satisfied or maxiter reached
        if isempty(ctr.ind)||iter==option.maxiter
            if iter==option.maxiter
                exitflag=0;status={'Maximum iterations reached'};
            else
                exitflag=1;status={'OK'};
            end
            output.iter=iter;output.status=status;output.drop=drop;
            return;
        end
    end
    %%% STEP2: find primal and dual step directions
    k=ctr.ind;%current constraint
    d=(Ai(k,:)*J)';%compute d
    List=iq:n;z = - J(:,List)*d(List);%update z  
    List=(1:iq-1);
    if ~isempty(List), r=-linsolve(R(List,List),d(List),opts);else, r=[];end%update r
    if norm(z) < option.tol
        TauPrimal = inf;
    else
        TauPrimal = (ctr.b - ctr.a'*x)/(ctr.a'*z);
    end 
    negative_dual_step_ind = find(r < -option.tol);
    if isempty(W) || isempty(negative_dual_step_ind)
        TauDual = inf;
        DualBlockingIndex = 0;
    else
        [TauDual,ind] = min(-u(W(negative_dual_step_ind))./r(negative_dual_step_ind));
        DualBlockingIndex = negative_dual_step_ind(ind);
    end
    %%% if TauPrimal == TauDual, TauPrimal is chosen
    [tau, TauCase] = min([TauPrimal, TauDual]);
    %%% determine new S-pair and make a step
    if TauPrimal == inf && TauDual == inf% no step in primal and dual space
        status ={'Infeasible Inequality Constraints'};
        exitflag=-2;
        output.iter=iter;output.status=status;output.drop=drop;
        return
    elseif TauPrimal == inf % step in dual space
        if ~isempty(W),u(W) = u(W) + tau * r;end
        ctr.u = ctr.u + tau;
        %%% drop blocking constraint
        [R,J,iq]=Downdate(R,J,iq,DualBlockingIndex);
        W(DualBlockingIndex) = [];%#ok      
        drop=drop+1;
        FullStep= false;  %%% goto STEP2
    else % step in primal and dual space  
        %TauPrimal < inf & TauDual < inf
        x = x + tau*z;
        if ~isempty(W),u(W) = u(W) + tau * r;end
        ctr.u = ctr.u + tau;
        if TauCase == 1 % TauPrimal: take a full step
            %%% TauPrimal <= TauDual
            %%% add selected primal constraint
            W = [W,ctr.ind];%#ok
            [R,J,iq]=UpdateRJ(R,J,d,iq,n);
            u(ctr.ind) = ctr.u;
            FullStep= true;  %%% goto STEP1       
        else % partial step
            %%% TauPrimal > TauDual & TauDual   <  inf
            %%% drop blocking constraint
            [R,J,iq]=Downdate(R,J,iq,DualBlockingIndex);
            W(DualBlockingIndex) = [];%#ok
            drop=drop+1;
            FullStep= false;  %%% goto STEP2
        end
    end
    iter=iter+1;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ctr=SelectConstraint(x,Ai,bi,nact,tol)
s=Ai*x-bi;s=s(nact);[val,ind]=max(s);
if val> tol, ind=nact(ind);else, ind=[];end
ctr.ind = ind;ctr.u   = 0;ctr.a  = Ai(ind,:)';ctr.b  = bi(ind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R,J,iq]=UpdateRJ(R,J,d,iq,n)
iqp1=iq+1;
[R(:,iq),u,b]=CHouseholder(d,iqp1:n,iq);
J=ApplyHouseholder(J,u,b,'R');
iq=iqp1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R,J,iq]=Downdate(R,J,iq,pos)
R(:,pos)=[];%delete column pos
%reduce the Hessenberg part of R to triangular
for k=pos:iq-2
    L=(k+1:iq-2);
    %Compute the Givens rotation such that R(k+1,k)=0
    [c,s,R(k,k),R(k+1,k)] = CGivens(R(k,k),R(k+1,k));%
    %Apply the transformation to the following columns from the left
    [R(k,L),R(k+1,L)] = ApplyGivens(c,s,R(k,L),R(k+1,L),'L');
    %apply the transposed transformation to J from the right
    [J(:,k),J(:,k+1)]=ApplyGivens(c,s,J(:,k),J(:,k+1),'RT');
end
iq=iq-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%