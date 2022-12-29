function [v1,v2]=ApplyGivens(c,s,v1,v2,option)
%[v1,v2]=ApplyGivens(c,s,v1,v2)
% Apply a Givens rotation Q=[c,-s;s,c] to vectors v1 and v2
% when option = 'L' ApplyGivens computes [v1;v2]=Q*[v1;v2]
%                                In this case v1 and v2 must be row vectors
% when option = 'R' ApplyGivens computes [v1,v2)=[v1,v2]*Q
%                                In this case v1 and v2 must be column vectors
%
% option = 'RT' or option = 'LT' do the same with Q=[c,s;-s,c] (opposite rotation)
% default value option='L'
%
% Author : Alain Barraud, copyright abc.consultant@wanadoo.fr 2003-2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<5, option='L';else, option=upper(option);end
if length(option)==2, s=-s;end
Q=[c,-s;s,c];
switch option(1)
    case 'L'%compute Q*[v1;v2],
        tmp=Q*[v1;v2];v1=tmp(1,:);v2=tmp(2,:);
    case 'R'%compute [v1,v2]*Q, 
        tmp=[v1,v2]*Q;v1=tmp(:,1);v2=tmp(:,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%