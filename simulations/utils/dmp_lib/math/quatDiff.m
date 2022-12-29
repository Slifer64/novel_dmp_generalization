function Qdiff = quatDiff(Q1, Q2)
%  quatDiff Calculate the product of two quaternions.
%   Qdiff = quatDiff(Q1, Q2) calculates the quaternion difference, Qdiff, 
%   for two given quaternions, Q1 and Q2. 
%   Qdiff represents the quaternion by which Q2 must be rotated to get Q1
%   in unit time, i.e. Q1 = Qdiff*Q2, where * is the quaternion product.
%   Each element of Q1 and Q2 must be a real number.  Additionally, Q1 and
%   Q2 have their scalar number as the first column.

Qdiff = quatProd(Q1, quatInv(Q2));

end