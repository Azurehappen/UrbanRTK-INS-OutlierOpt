function C_e_e = computeRotForEcefToEnu(lla_deg)
% Compute the Earth to Geographic rotation matrix given lat long.     
%                                                                          
% Parameters:                                                              
%   llh_b:      3-by-1 vector, [Lat (rad), Long (rad), Height (m)] position.                           
%                                                                          
% Return values:                                                                                         
%   C_e_n:      3-by-3 rotation matrix from ECEF to ENU                                                     
%                                                                          
% Reference:                                                                                                                             
%   - Paul Groves "Principles of GNSS, Inertial, and Multisensor
%     Integrated Navigation Systems," Second Edition.

% calculations
%-------------------------------------------------------------------------%
lat_rad = deg2rad(lla_deg(1));
lng_rad = deg2rad(lla_deg(2));
sin_lat= sin(lat_rad);
cos_lat= cos(lat_rad);
sin_lng= sin(lng_rad);
cos_lng= cos(lng_rad);
C_e_e = [-sin_lng          cos_lng           0;
         -sin_lat*cos_lng -sin_lat*sin_lng   cos_lat;
         cos_lat*cos_lng cos_lat*sin_lng  sin_lat];