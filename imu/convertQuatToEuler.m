function [roll_phi_rad,pitch_theta_rad,yaw_psi_rad] = convertQuatToEuler(q_e_b)
% Parameters:                                                              
%   q_a_b:	4-by-1 quaternion vector = [a bi cj dk]                                    
%                                                                          
% Return values:                                                           
%   R_a_b:	3-by-3 rotation matrix describing transformation from Body to 
%           Nav. This is actually valid for any alpha to beta rotation. 

% Convert quaternion to Rotation matrix representing Body to Nav, i.e. R_b_n
% Farrell D.2.2, works for any normal quaternion
% roll: phi, pitch: theta, yaw: psi
if norm(q_e_b)~=0.0
    q_e_b = q_e_b/norm(q_e_b);
    b1 = q_e_b(1);
    b2 = q_e_b(2);
    b3 = q_e_b(3);
    b4 = q_e_b(4);
    % The euler angle is from body to ECEF
    sin_theta = -2*(b2*b4+b1*b3);
    pitch_theta_rad = asin(sin_theta);
    roll_phi_rad = limit(atan2(2*(b3*b4-b1*b2), 1-2*(b2^2+b3^2)));
    yaw_psi_rad = limit(atan2(2*(b2*b3-b1*b4), 1-2*(b3^2+b4^2)));
else
    R_a_b = eye(3); % fault condition
    error('Norm b=0 in convertQuatToRot()');
end

function [y]=limit(x)

while x>pi
    x=x-2*pi;
end
while x<=-pi
    x=x+2*pi;
end
y=x;