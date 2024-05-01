function rot_b2n = computeDirectionCos(phi,theta,psi)

% Define plane rotation, Eqns. 2.24-2.26 in Farrell Book
rot_x = [1,0,0;
    0,cos(phi),sin(phi);
    0,-sin(phi),cos(phi)]; % roll
rot_y = [cos(theta),0,-sin(theta);
    0,1,0;
    sin(theta),0,cos(theta)]; % pitch
rot_z = [cos(psi),sin(psi),0;
    -sin(psi),cos(psi),0;
    0,0,1]; % yaw

rot_b2n = rot_z'*rot_y'*rot_x'; % Section 10.3

