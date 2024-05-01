function [elev, az] = sat_elev_azimuth(p,u_ecef, s_ecef) 
% Compute the azimuth angle and elevation angle between user and satellite
% Input:
%       u_ecef: 3 by 1 vector, user position in ECEF
%       s_ecef: 3 by 1 vector, satellite position in ECEF

%geodetic position (WGS-84)
[~, ~, ~, R_e2l, ~, ~] = ecef2llh(p,u_ecef);
% user and satellite positions in local level coordinate system
r_u_level = R_e2l*u_ecef;
r_s_level = R_e2l*s_ecef;
dr = r_s_level - r_u_level;

% elevation angle between user    and satellite (semi-circles)
elev = atan2(-dr(3), norm(dr(1:2)));

% azimuth angle between user and satellite, CW positive from true north (semi-circles)
az = atan2(dr(2), dr(1));
end

