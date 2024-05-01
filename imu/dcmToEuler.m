function [roll_phi_rad,pitch_theta_rad,yaw_psi_rad] = dcmToEuler(R_b2n)
% Convert DCM to roll,pitch,yaw representing Body to Nav, i.e. R_b_n
% Farrell 2.45, works for any normal quaternion
% roll: phi, pitch: theta, yaw: psi

pitch_theta_rad = -atan(R_b2n(3,1)/sqrt(1-(R_b2n(3,1))^2));
roll_phi_rad = limit(atan2(R_b2n(3,2), R_b2n(3,3)));
yaw_psi_rad = limit(atan2(R_b2n(2,1), R_b2n(1,1)));

function [y]=limit(x)

while x>pi
    x=x-2*pi;
end
while x<=-pi
    x=x+2*pi;
end
y=x;