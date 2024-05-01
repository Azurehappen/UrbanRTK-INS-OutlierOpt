function [R_a2b] = convertQuatToRot(q_a2b)
% Parameters:                                                              
%   q_a_b:	4-by-1 quaternion vector = [a bi cj dk]                                    
%                                                                          
% Return values:                                                           
%   R_a_b:	3-by-3 rotation matrix describing transformation from Body to 
%           Nav. This is actually valid for any alpha to beta rotation. 

% Convert quaternion to Rotation matrix representing Body to Nav, i.e. R_b_n
% Farrell eqn. D13, works for any normal quaternion
if norm(q_a2b)~=0.0
    q_a2b = q_a2b/norm(q_a2b);
    b1 = q_a2b(1);
    bv = q_a2b(2:4);
    R_a2b = (b1^2 - bv'*bv)*eye(3,3) + 2*(bv*bv') + 2*b1*vectorSkewSymMat(bv); % Farrell eq. D.14
    % Note: parenthesize Bv*Bv' to ensure result is Hermitian 
else
    R_a2b = eye(3); % fault condition
    error('Norm b=0 in convertQuatToRot()');
end