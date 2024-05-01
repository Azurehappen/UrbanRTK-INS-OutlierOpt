function [q] = quatMult(a,b)                          
% Multiply two quaternions.                                                
%          
% Syntax:
%   [q] = quatMult(a,b)                                                
%                                                                          
% Parameters:                                                              
%   a:      4-by-1 quaternion vector 1                                     
%   b:      4-by-1 quaternion vector 2                                     
%                                                                          
% Return values:                                                           
%   q:      4-by-1 quaternion vector multiplied                            
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
% ID:                   $Id: quatMult.m 391 2017-06-30 22:13:40Z proysdon $
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
% Farrell eqn. D.3
A = [a(1) -a(2) -a(3) -a(4)
     a(2)  a(1) -a(4)  a(3)
     a(3)  a(4)  a(1) -a(2)
     a(4) -a(3)  a(2)  a(1)];
q = A*b;

% this is Madgwick's equation
% q(1,1) = a(1).*b(1)-a(2).*b(2)-a(3).*b(3)-a(4).*b(4);
% q(2,1) = a(1).*b(2)+a(2).*b(1)+a(3).*b(4)-a(4).*b(3);
% q(3,1) = a(1).*b(3)-a(2).*b(4)+a(3).*b(1)+a(4).*b(2);
% q(4,1) = a(1).*b(4)+a(2).*b(3)-a(3).*b(2)+a(4).*b(1);

end