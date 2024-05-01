function [llh_b, R_e_n, g_b, omega_ie_e] = ecef2llh_iter(r_eb_e)
% Convert ECEF to LLH representation, using iteration.                                      
%                                                                          
% Syntax:                                                                  
%   [llh_b, R_e_n, g_b, omega_ie_e] = ecef2llh(r_eb_e)            
%                                                                          
% Parameters:                                                              
%   r_eb_e:     3-by-1 vector, [x,y,z] position ECEF (m)                            
%                                                                          
% Return values:                                                           
%   llh_b:      3-by-1 vector, [Lat (rad), Long (rad), Height (m)] position.                                  
%   R_e_n:      3-by-3 rotation matrix from ECEF to NED                    
%   g_l:        3-by-1 vector, local gravity (m/s/s)                                      
%   omega_ie_e:	3-by-1 vector, Earth angular rate (rad/s)                                       
%                                                                          
% Reference:                                                               
%   - IS-GPS-200 20.3.3.4.3.3.2 Earth-centered, inertial (ECI) coordinate    
%     system                                                                 
%   - Paul Groves "Principles of GNSS, Inertial, and Multisensor
%     Integrated Navigation Systems," Second Edition.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author(s):            P.F. Roysdon                                       
% Last committed:       $Revision: 439 $                                   
% Last changed by:      $Author: proysdon $                                  
% Last changed date:    $Date: 2017-08-22 15:01:39 -0700 (Tue, 22 Aug 2017) $
% ID:                   $Id: ecef2llh_iter.m 439 2017-08-22 22:01:39Z proysdon $                                                  
% email:                software@aidednav.com
% Website:              http://www.aidednav.com
% 
% Copyright ©, 2016 Aided Nav,  All rights reserved.
%
% Approved for use by University of California Riverside, Department of 
% Electrical Engineering, Control and Robotics Lab. 
%                                                                          
% This program carries no warranty, not even the implied                   
% warranty of merchantability or fitness for a particular purpose.         
% 
% Please email bug reports or suggestions for improvements to:
% software@aidednav.com or proysdon@ece.ucr.edu
%


% calculations
%-------------------------------------------------------------------------%
% WGS84 Earth model
g = 9.80665;                    % standard gravity
omega_ie = 7.29211e-5;          % angular rate
a = 6378137.0;                  % semi-major axis length, meters
b = 6356752.31424518;           % semi-minor axis length, meters
f = 1.0/298.257223563;          % flatness of ellipsoid
e = sqrt((a*a-b*b)/(a*a));      % first eccentricity of ellipsoid
ep = sqrt((a*a-b*b)/(b*b));     % second eccentricity of ellipsoid

% convert xyz to llh
h_b = 0;
N = a;
p = norm(r_eb_e(1:2,1));
z = r_eb_e(3,1);
L_b = 0;

for iter = 0:10
	s_L_b = saturate(z / (N*(1.0 - e^2) + h_b), 1);
	L_b = atan2((z + e^2*N*s_L_b), p);
	N = a / sqrt(1.0 - (e*s_L_b)^2);
	h_b = p / cos(L_b) - N;
end

lambda_b = atan2(r_eb_e(2,1), r_eb_e(1,1));

% LLH
llh_b      = zeros(3,1);
llh_b(1,1) = L_b; % latitude (rad)
llh_b(2,1) = lambda_b; % longitude (rad)
llh_b(3,1) = h_b; % height (m)

% sin|cos conversions (do this once)
sin_lat = sin(L_b);
cos_lat = cos(L_b);
sin_lon = sin(lambda_b);
cos_lon = cos(lambda_b);

% earth-to-local-level (nav) rotation.  Note: local-level also  called NED
R_e_n = [-sin_lat*cos_lon, -sin_lat*sin_lon,  cos_lat;
		 -sin_lon,          cos_lon,          0;
		 -cos_lat*cos_lon, -cos_lat*sin_lon, -sin_lat];

% local gravity vector
g_b    = zeros(3,1);
g_b(3) = 9.780318 * (1.0 + 0.0053024*sin_lat^2 ...
          - 0.0000058*sin(2*L_b)^2) - 3.086e-6*h_b;

% Earth angular rate vector
omega_ie_e = [0; 0; omega_ie];

