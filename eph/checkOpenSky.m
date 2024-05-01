function is_open_sky = checkOpenSky(gps_range, gps_sat_pos, user_pos)
% Check the open sky condition based on SAE J2945/1_202004
% 1. No view obstructions external to the vehicle can be seen, from the
% point of view of the GNSS antennas of reference device and device under
% test, starting from 5 degrees above the horizontal plane (the elevation
% mask) containing the antenna phase center, in all directions around
% the antenna.
% 2. The number of healthy satellites used, as reported by the reference
% device, for GPS satellites only, is greater than or equal to seven.
% 3. The HDOP, as reported by the reference device, for GPS satellites,
% is less than or equal to 1.5, and VDOP is less than or equal to three. 
%
% Input: gps_range: GPS code measurement with sat above 5 deg elevation.
%        gps_sat_pos: GPS sat posion in Er frame (Rotated).
%        user_pos: Coarse user position, a 3x1 vector.


if length(gps_range) < 7
    is_open_sky = false;
    return;
end
num_y = length(gps_range);
range = zeros(num_y,1);
H = zeros(num_y,3+1);
H(:,4) = 1;
for j=1:num_y
    range(j)=norm(gps_sat_pos(:,j)-user_pos);
    % Line of sight
    H(j,1:3) = (user_pos-gps_sat_pos(:,j))'/range(j);
end
lla_deg = ecef2lla(user_pos', 'WGS84');
R_e2g=computeRotForEcefToNed(lla_deg);
R = [R_e2g, zeros(3,1);
    zeros(1,3),1];
Q = (H'*H)^(-1);
Q_ned = R*Q*R';
HDOP = sqrt(Q_ned(1,1)+Q_ned(2,2));
VDOP = sqrt(Q_ned(3,3));
if HDOP <= 1.5 && VDOP <= 3
    is_open_sky = true;
    return;
end
is_open_sky = false;
end