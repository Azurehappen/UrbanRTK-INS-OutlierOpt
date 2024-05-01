function [lat_f, lon_f]= pierce_point_calc(lat_0,lon_0,azim,elev,p)
%%%% input:
%%% lat_0: user latitude position
%%% lon_0: user longitude position
%%%% azim: azimuth angle in radians
%%%% elev: elevation angle in radians
%%% r: user radius
%%% p: parameter constant

%%%% output:
%%%%lat_f: latitude of pierce point
%%%%lon_f: longitude of pierce point

%%%%% reference:
%%%%  COMPARATIVE STUDY OF METHODS FOR CALCULATING IONOSPHERIC POINTS AND
%%%%  DESCRIBING THE GNSS SIGNAL PATH (Prol-2017)

%%%% initialization
lat_r=lat_0;
lon_r=lon_0;
elev_r=elev;
azim_r=azim;

%%% From (Prol-2017)
%%% eqn 3
psi_r=pi/2-elev-asin(p.Re*cos(elev)/(p.Re+p.h_iono));
%%% eqn 1
rr1=sin(lat_r)*cos(psi_r)+cos(lat_r)*sin(psi_r)*cos(azim_r);

%%% latitude of the ionosphere piercing point
phi_ip=asin(rr1);

% %%%% longitude of the ionosphere piercing point
% if ((phi_ip >  70*pi/180) & (tan(psi_r)*cos(azim)      > tan((pi/2) - lat_r))) | ...
%    ((phi_ip< -70*pi/180) & (tan(psi_r)*cos(azim + pi) > tan((pi/2) + lat_r)))
%       lambda_ip = lon_r + pi - asin((sin(psi_r)*sin(azim))/(cos(phi_ip)));
% else
%     lambda_ip = lon_r + asin((sin(psi_r)*sin(azim))/(cos(phi_ip)));
% end

%%% eqn 2
lambda_ip=lon_r+asin(sin(psi_r)*sin(azim_r)/cos(phi_ip));

%%% convert it from radians to degree
lat_f=phi_ip*180/pi;
lon_f=lambda_ip*180/pi;




