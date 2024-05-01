function [vel,clock_rate] = velSolver(pos,cpt)
 
yv = cpt.doppler;
num = length(yv); % The number of measurement
H = zeros(num,4);
R = zeros(num,1); 
s_pos_ecef = cpt.sat_pos_Rcorr;
s_v_ecef = cpt.sat_v_Rcorr;
res_v = zeros(length(yv),1);
rate = zeros(length(yv),1);
for j=1:num
    R(j)=norm(s_pos_ecef(:,j)-pos);
    los= (pos-s_pos_ecef(:,j))'/R(j);
    H(j,:)=[los 1];
    rate(j) = los*s_v_ecef(:,j);
    res_v(j) = yv(j) + los*s_v_ecef(:,j);
end
vel = (H'*H)^(-1)*H'*res_v;

clock_rate = vel(4);
vel = vel(1:3);
% wgs84 = wgs84Ellipsoid('meter');
% [xEast,yNorth,zUp] = ecef2enu(vel(1),vel(2),vel(3),lat0,lon0,h0,wgs84)
end