function [c,s,aa,bb] = CGivens(a,b)
%[c,s,aa,bb]=CGivens(a,b)
% Compute a Givens rotation Q=[c,-s;s,c] such that the transformed 
% vector [aa;bb] = Q[a;b] = [w 0] with abs(w)=norm([a;b])
%
% Author : Alain Barraud, copyright abc.consultant@wanadoo.fr 2003-2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if b~= 0
    aa = hypot(a,b);bb=0;
    c = a/aa;s = -b/aa;
else
    aa=a;bb=0;c = 1;s = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%