function [t_out] = limit_tgps(t)                                  
% Limit GPS time (in seconds) to 1-week according to GPS-ICD-200.          
%                                                                          
% Syntax:                                                                  
%   [t] = limit_tgps(t)                                                    
%                                                                          
% Parameters:                                                              
%   t:      time of week                                                   
%                                                                          
% Return values:                                                           
%   t_out:  time of week                                                   
%                                                                          


% calculations
%-------------------------------------------------------------------------%
excd = find(t > 302400);
less = find(t < -302400);
if ~isempty(excd)
    t(excd) = t(excd) - 604800;
elseif ~isempty(less)
    t(less) = t(less) + 604800;
end
t_out = t;

end