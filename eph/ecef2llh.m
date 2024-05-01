function [lat, lon, height, R_e2l, gP_l, omega_ie_e] = ecef2llh(p,xyz) 
% Convert ECEF to LLH representation.                                      
%                                                                          
% Syntax:                                                                  
%   [lat, lon, height, R_e2l, gP_l, omega_ie_e] = ecef2llh(xyz)            
%                                                                          
% Parameters:                                                              
%   xyz:        3-by-1 position vector (m) ECEF                            
%                                                                          
% Return values:                                                           
%   lat:        position in Latitude (radians)                                 
%   lon:        position in Longitude (radians)                                
%   height:     Altitude above user datum ellipsoid  (m)                     
%   R_e2l:      3-by-3 rotation matrix from ECEF to ECI                    
%   gP_l:       local gravity (m/s/s)                                      
%   omega_ie_e:	angular rate (rad/s)                                       
%                                                                          
% Reference:                                                               
%   IS-GPS-200 20.3.3.4.3.3.2 Earth-centered, inertial (ECI) coordinate    
%   system                                                                                                                                    


% computations
%--------------------------------------------------------------%
omega_ie = p.omge;
a = p.a;
e = p.e;
% convert xyz to llh
h = 0;
N = a;
pp = norm(xyz(1:2));
z = xyz(3);
lambda = 0;
limit = 1;
for iter = 0:10
    % Saturate signal
    s_lambda = z/(N*(1.0 - e^2) + h);
    if s_lambda > limit
        s_lambda = limit;
    elseif s_lambda < -limit
        s_lambda = -limit;
    end

	lambda = atan2((z + e^2*N*s_lambda), pp);
	N = a / sqrt(1.0 - (e*s_lambda)^2);
	h = pp / cos(lambda) - N;
end

lat = lambda;
lon = atan2(xyz(2), xyz(1));
height = h;
% sin|cos conversions (do this once)
slat = sin(lat);
clat = cos(lat);
slon = sin(lon);
clon = cos(lon);

% earth-to-level rotation
R_e2l = [-slat*clon, -slat*slon,  clat;
		 -slon,       clon,       0;
		 -clat*clon, -clat*slon, -slat];

gP_l = zeros(3,1);
gP_l(3) = 9.780318 * (1.0 + 0.0053024*slat^2 - 0.0000058*sin(2*lat)^2) - 3.086e-6*height;

omega_ie_e = [0; 0; omega_ie];


end