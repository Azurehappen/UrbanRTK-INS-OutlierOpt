function [qc] = quatConj(q)
% Compute the quaternion conjugate
% 
% [qc] = quatConj(q)
% 
% Parameters:
%   q:      4-by-1 quaternion vector
%   
% Return values:
% 	qc:     4-by-1 quaternion vector conjugate
% 
% Reference:
%	J. A. Farrell, "Aided Navigation: GPS with High Rate Sensors", 
%   McGraw Hill, 2008.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author(s):            P.F. Roysdon                                       
% Last committed:       $Revision: 432 $                                   
% Last changed by:      $Author: proysdon $                                  
% Last changed date:    $Date: 2017-08-10 10:00:29 -0700 (Thu, 10 Aug 2017) $
% ID:                   $Id: quatConj.m 432 2017-08-10 17:00:29Z proysdon $                                                  
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
% Farrell eqn. D.2
qc = [q(1); -q(2:4)];

