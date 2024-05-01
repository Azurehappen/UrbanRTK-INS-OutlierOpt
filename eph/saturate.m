function [x]=saturate(x, limit)                                   
% Saturate signal.                                                         
%                                                                          
% Syntax:                                                                  
%   [x]=saturate(x, limit)                                                 
%                                                                          
% Parameters:                                                              
%   x:      signal to saturate                                             
%   limit:  saturation limit                                               
%                                                                          
% Return values:                                                           
%   x:      saturated signal                                               
%                                                                          
% Reference:                                                               
%   IS-GPS-200                                                             
%
% Other m-files required: hmat
% Subfunctions: none
% MAT-files required: none
%

% Author(s):            P.F. Roysdon                                       
% Last committed:       $Revision: 22 $                                   
% Last changed by:      $Author: proysdon $                                  
% Last changed date:    $Date: 2016-08-07 19:05:13 -0700 (Sun, 07 Aug 2016) $
% ID:                   $Id: saturate.m 22 2016-08-08 02:05:13Z proysdon $                                                  
% email:                software@aidednav.com
% Website:              http://www.aidednav.com
% 
% Copyright ©, 2016 Aided Nav,  All rights reserved.
%
% Approved for use by University of California Riverside, Department of 
% Electrical Engineering, Robotics and Control Lab. 
%                                                                          
% This program carries no warranty, not even the implied                   
% warranty of merchantability or fitness for a particular purpose.         
% 
% Please email bug reports or suggestions for improvements to:
% software@aidednav.com or proysdon@ee.ucr.edu
%                                                                     


% computations
%--------------------------------------------------------------%
i = find(x > limit);
if ~isempty(i),
    x(i) = limit;
end
i = find(x < -limit);
if ~isempty(i),
    x(i) = -limit;
end


end