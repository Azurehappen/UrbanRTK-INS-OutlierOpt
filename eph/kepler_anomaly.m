function [E_k] = kepler_anomaly(ecc, M_k)                         
% Compute the GPS Kepler anomaly.                                          
%                                                                          
% Description:                                                             
% This function solves Keplers equation for eccentric anaomaly E in terms  
% of mean anomaly M and eccentricity ecc for low eccentricity elliptical   
% orbits                                                                   
%                                                                          
% Syntax:                                                                  
%   [E_k] = kepler_anomaly(ecc, M_k)                                       
%                                                                          
% Parameters:                                                              
%   ecc:    orbit eccentricity  (unitless)                                 
%   M_k:    mean anomaly in radians                                        
%                                                                          
% Return values:                                                           
%	E_k:    Eccentric anomaly                                              
%                                                                          
% Reference:                                                               
%   IS-GPS-200                                                             
%
% Other m-files required: hmat
% Subfunctions: none
% MAT-files required: none
%

% Author(s):            P.F. Roysdon
% Last committed:       $Revision: 391 $
% Last changed by:      $Author: proysdon $
% Last changed date:    $Date: 2017-06-30 15:13:40 -0700 (Fri, 30 Jun 2017) $
% ID:                   $Id: kepler_anomaly.m 391 2017-06-30 22:13:40Z proysdon $
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
% set tolerance
tol = 1e-13;
max_iter = 30;

% initial guess
%E_k = M_k + ecc*sin(M_k) / (1 - sin(M_k + ecc) + sin(M_k));
E_k = M_k;
% improve guess with Newton-Raphson Iteration
for iter = 1:max_iter
	Delta_E = (E_k - ecc*sin(E_k) - M_k) / (1 - ecc*cos(E_k));

    absDelta_E = abs(Delta_E);
    
    E_k = E_k - Delta_E;
    
	if all(absDelta_E < tol)
		break;
	end
	
    % if (i == max_iter)
    %     fprintf(1,'WARNING: Kepler anomaly does not converge.\n');
    % end
end


end