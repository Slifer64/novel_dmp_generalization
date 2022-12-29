function Qconj = quatConj(Q)
%  quatInv Calculate the inverse of a quaternion.
%   invQ = quatInv(Q) calculates the quaternion inverse, invQ, of Q. 
%   Q must have the form Q = [n e] where n is the scalar and e the vector
%   part.
%   Note: Q does not need to be a unit quaternion.

Qconj = zeros(size(Q));

Qconj(1) = Q(1);
Qconj(2:4) = -Q(2:4);

end