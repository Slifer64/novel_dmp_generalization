function qmat = quat2mat(Q)
%  quatProd Calculate the product of two quaternions.
%   Q12 = quatProd(Q1, Q2) calculates the quaternion product, Q12, for two
%   given quaternions, Q1 and Q2. Each element of Q1 and Q2 must be a
%   real number.  Additionally, Q1 and Q2 have their scalar number as the first 
%   column.
%   Note: Quaternion multiplication is not commutative.

qmat = zeros(4,4);

w = Q(1);
v = Q(2:4);

qmat(:,1) = Q;
qmat(1,2:4) = -v;
qmat(2:4,2:4) = w*eye(3,3) + vector2ssMatrix(v);

end

function ssMat = vector2ssMatrix(p)

ssMat(1,1) = 0;
ssMat(2,2) = 0;
ssMat(3,3) = 0;
ssMat(1,2) = -p(3);
ssMat(2,1) = p(3);
ssMat(1,3) = p(2);
ssMat(3,1) = -p(2);
ssMat(3,2) = p(1);
ssMat(2,3) = -p(1);

end