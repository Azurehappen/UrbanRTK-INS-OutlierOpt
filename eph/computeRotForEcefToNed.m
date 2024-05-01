function [C_e_n] = computeRotForEcefToNed(lla_deg)
% Compute the Earth to Geographic rotation matrix given lat long.           
%                                                                          
% Parameters:                                                              
%   llh_b:      3-by-1 vector, [Lat (deg), Long (deg), Height (m)] position.                           
%                                                                          
% Return values:                                                                                         
%   C_e_n:      3-by-3 rotation matrix from ECEF to NED                                                         
%                                                                          
% Reference:                                                                                                                             
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
% ID:                   $Id: ll2R.m 439 2017-08-22 22:01:39Z proysdon $                                                  
% email:                software@aidednav.com
% Website:              http://www.aidednav.com
% 
% Copyright ?, 2016 Aided Nav,  All rights reserved.
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
lat_rad = deg2rad(lla_deg(1));
lng_rad = deg2rad(lla_deg(2));
sin_lat= sin(lat_rad);
cos_lat= cos(lat_rad);
sin_lng= sin(lng_rad);
cos_lng= cos(lng_rad);
C_e_n = [   -sin_lat*cos_lng -sin_lat*sin_lng   cos_lat; % Farrell eqn. 2.32
            -sin_lng          cos_lng           0;
            -cos_lat*cos_lng -cos_lat*sin_lng  -sin_lat];