function quat = directionCosToQuat(rot_b2e)
% Convert a rotation matrix to a normalized quaternion

% Eqn. D.15 in Farrell Book
quat = zeros(4,1);
quat(1) = 0.5*sqrt(1+rot_b2e(1,1)+rot_b2e(2,2)+rot_b2e(3,3));
quat(2) = (rot_b2e(3,2)-rot_b2e(2,3))/(4*quat(1));
quat(2) = (rot_b2e(1,3)-rot_b2e(3,1))/(4*quat(1));
quat(2) = (rot_b2e(2,1)-rot_b2e(1,2))/(4*quat(1));
quat = quat/norm(quat);
end