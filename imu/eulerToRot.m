function [R_e_b] = eulerToRot(roll_phi,pitch_theta,yaw_psi)                 
% Convert Euler angles to a rotation matrix, E = [roll, pitch, yaw]                                                                                     
%        
% Syntax:
%   [R_n_b] = euler2R(E_b)                                                 
%                                                                          
% Parameters:                                                              
%   E_b:    3-by-1 vector, [roll, pitch, yaw] (rad)                             
%                                                                          
% Return values:                                                           
%   R_n_b:  3-by-3 rotation matrix describing transformation from beta to alpha                                                                                                             
% 
% Reference:
%	Groves, P. "Principles of GNSS, Inertial, and Multisensor
%	Integrated Navigation Systems," Second Edition.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author(s):            P.F. Roysdon                                       
% Last committed:       $Revision: 432 $                                   
% Last changed by:      $Author: proysdon $                                  
% Last changed date:    $Date: 2017-08-10 10:00:29 -0700 (Thu, 10 Aug 2017) $
% ID:                   $Id: euler2R.m 432 2017-08-10 17:00:29Z proysdon $                                                  
% email:                software@aidednav.com
% Website:              http://www.aidednav.com
% 
% Copyright �, 2016 Aided Nav,  All rights reserved.
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
% cPhi = cos(E_b(1)); % phi = roll
% sPhi = sin(E_b(1));
% cThe = cos(E_b(2)); % theta = pitch
% sThe = sin(E_b(2));
% cPsi = cos(E_b(3)); % psi = yaw
% sPsi = sin(E_b(3));

% Rotation matrix representing Nav to Body, i.e. C_n_b
% R_n_b =[ cPsi*cThe                   sPsi*cThe                   -sThe
%         (-sPsi*cPhi+cPsi*sThe*sPhi) ( cPsi*cPhi+sPsi*sThe*sPhi)  (cThe*sPhi)
%         ( sPsi*sPhi+cPsi*sThe*cPhi) (-cPsi*sPhi+sPsi*sThe*cPhi)  (cThe*cPhi)];  % Farrell eqn 2.43

cpsi = cos(yaw_psi); 
spsi = sin(yaw_psi);
cthe = cos(pitch_theta); 
sthe = sin(pitch_theta);
cphi = cos(roll_phi); 
sphi = sin(roll_phi);

C1 = [ cpsi  spsi  0   ; ...
      -spsi  cpsi  0   ; ...
       0     0     1   ];
C2 = [ cthe  0    -sthe; ...
       0     1     0   ; ...
       sthe  0     cthe];
C3 = [ 1     0     0   ; ...
       0     cphi  sphi; ...
       0    -sphi  cphi]; 

R_e_b = C3*C2*C1;
