function iono_delay = klobuchar_model(p,iono,elev,az,t_gps)
% Compute GPS satellite ionosphere corrections according to IS-GPS-200.    
%                                                                          
% Description:                                                             
% This function computes the single frequency ionospheric delay correction 
% in seconds as specified in Figure 20-4 in ICD-GPS-200                    
%                                                                          
% Syntax:                                                                  
%   [dt_iono, az, elev, eFlg] = sat_iono(p, eph, iono, t_gps, r_u_ecefeci, r_s_eci)
%                                                                          
% Parameters:                                                              
%   p:           structure containing parameters (settings)                      
%   iono:        structure containing satellite ionosphere values          
%                 - alpha0, alpha1, alpha2, alpha3 iono model parameters   
%                   from 50 bps message                                    
%                 - beta0, beta1, beta2, beta3 iono model parameters       
%                   from 50 bps message
%   az:         satellite azimuth (rad)                                    
%   elev:       satellite elevation (rad)
%   t_gps:       GPS time of signal reception (s)                                                                                                  
% Return values:                                                           
%   iono_delay:    iono delay (m)                                                                                     
%                                                                          
% Reference:                                                               
%   - IS-GPS-200                                                           
%   - J. A. Farrell, "Aided Navigation: GPS with High Rate Sensors",       
%     McGraw Hill, 2008.
dt_iono = 0;

% geodetic position in semi-circles
phi_u = p.lat/p.pi;
lambda_u = p.lon/p.pi;

% earth's central angle between user position and the earth project of ionospheric intersection point (semi-circles)
Psi = 0.0137/((elev/p.pi) + 0.11) - 0.022;

% geodetic latitude of the earth projection of the ionospheric intersection point (semi-circles)
phi_i = saturate(phi_u + Psi*cos(az), 0.416);

% geodetic longitude of the earth projection of the ionospheric intersection point (semi-circles)
lambda_i = lambda_u + Psi*sin(az)/cos(phi_i*p.pi);

% geomagnetic latitude of the earth projection of the ionospheric intersection point (semi-circles)
phi_m = phi_i + 0.064*cos((lambda_i - 1.617)*p.pi);

% local time (sec)
t = 4.32e4*lambda_i + t_gps;
if t >= 86400 
     t = t - 86400;
elseif t < -86400
     t = t + 86400;
end

% obliquity factor (dimensionless)
F = 1.0 + 16.0*(0.53 - (elev/p.pi))^3;

% period
PER = iono.ionoBeta(1)*phi_m^0 + iono.ionoBeta(2)*phi_m^1 ...
      + iono.ionoBeta(3)*phi_m^2 + iono.ionoBeta(4)*phi_m^3;
if PER < 72000
    PER = 72000;
end

% phase (radians)
x = 2*p.pi*(t - 50400)/PER;

% amplitude of delay
AMP = iono.ionoAlpha(1)*phi_m^0 + iono.ionoAlpha(2)*phi_m^1 ...
      + iono.ionoAlpha(3)*phi_m^2 + iono.ionoAlpha(4)*phi_m^3;
    if AMP < 0
        AMP = 0;
    end

% L1 ionosphere correction (s)
if abs(x) < 1.57
    dt_iono = F * (5e-9 + AMP*(1 - x^2/2 + x^4/24));
else
    dt_iono = F * (5e-9);
end
iono_delay = p.c * dt_iono;