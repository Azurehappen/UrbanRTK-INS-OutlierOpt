function [qinv] = quatInv(q)                                  
% Computes the quaternion inverse.
% 
%     
% Syntax:
%   [qinv] = quatinv(q)                                                   
%                                                                          
% Parameters:                                                              
%   q:      4-by-1 quaternion vector                                       
%                                                                          
% Return values:                                                           
%   qinv:   4-by-1 quaternion vector inverse                               
%                                                                          
% Reference:                                                               
%   J. A. Farrell, "Aided Navigation: GPS with High Rate Sensors",         
%   McGraw Hill, 2008.                                                     
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author(s):            P.F. Roysdon
% Last committed:       $Revision: 391 $
% Last changed by:      $Author: proysdon $
% Last changed date:    $Date: 2017-06-30 15:13:40 -0700 (Fri, 30 Jun 2017) $
% ID:                   $Id: quatInv.m 391 2017-06-30 22:13:40Z proysdon $
% email:                software@aidednav.com
% Website:              http://www.aidednav.com
% 
% Copyright, 2017 Aided Nav,  All rights reserved.
% 
% Approved for use by University of California Riverside, Department of 
% Electrical Engineering, Robotics and Control Lab. 
% 
% This program carries no warranty, not even the implied warranty of
% merchantability or fitness for a particular purpose.
% 
% Please email bug reports or suggestions for improvements to:
% software@aidednav.com
%                                                                     


% calculations
%-------------------------------------------------------------------------%
% Farrell eqn. D.5-6
qinv = quatConj(q)/norm(q);
