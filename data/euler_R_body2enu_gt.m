function r = euler_R_body2enu_gt(h_deg, p_deg, r_deg)
% Conversion info from Texus GT
heading = deg2rad(h_deg);
pitch = deg2rad(p_deg);
roll = deg2rad(r_deg);
R1 = [cos(heading),sin(heading),0;
    -sin(heading),cos(heading),0;
    0,0,1];
R2 = [1,0,0;
    0,cos(pitch),-sin(pitch);
    0,sin(pitch),cos(pitch)];
R3 = [cos(roll),0,sin(roll);
    0,1,0;
    -sin(roll),0,cos(roll)];
r = R1*R2*R3;
end