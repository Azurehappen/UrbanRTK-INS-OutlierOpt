function [q_e2b] = R2quat(R_e2b)
% Convert a rotation matrix to a normalized quaternion
%
% Syntax:
%   [q_b] = R2quat(R_b_n)
%
% Parameters:
%   R_b_n:  3-by-3 rotation matrix describing transformation from Body to 
%           Nav.  This is actually valid for any alpha to beta rotation. 
%
% Return values:
% 	q_b:    4-by-1 quaternion vector = [a bi cj dk]
%                                                                          
% Reference:  
% - J. A. Farrell, "Aided Navigation: GPS with High Rate Sensors", 
%   McGraw Hill, 2008.
% - Paul Groves "Principles of GNSS, Inertial, and Multisensor
%   Integrated Navigation Systems," Second Edition.
%


% calculations
%-------------------------------------------------------------------------%
b1 = 0.5*sqrt(1+R_e2b(1,1)+R_e2b(2,2)+R_e2b(3,3));
b2 = (R_e2b(3,2)-R_e2b(2,3))/(4*b1);
b3 = (R_e2b(1,3)-R_e2b(3,1))/(4*b1);
b4 = (R_e2b(2,1)-R_e2b(1,2))/(4*b1);
q_e2b = [b1;b2;b3;b4];
q_e2b = q_e2b/norm(q_e2b);

