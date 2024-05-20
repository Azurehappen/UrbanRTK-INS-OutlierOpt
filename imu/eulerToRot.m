function [R_n_b] = eulerToRot(roll_phi,pitch_theta,yaw_psi)                 
% Convert Euler angles to a rotation matrix, E = [roll, pitch, yaw]                                                                                     
%        
% Syntax:
%   [R_n_b] = eulerToRot(E_b)                                                 
%                                                                          
% Parameters:                                                              
%   E_b:    3-by-1 vector, [roll, pitch, yaw] (rad)                             
%                                                                          
% Return values:                                                           
%   R_n_b:  3-by-3 rotation matrix describing transformation from beta to alpha                                                                                                             
% 
% Farrell eqn 2.43

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

R_n_b = C3*C2*C1;
