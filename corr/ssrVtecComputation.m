function [iono_corr, mapping_m] = ssrVtecComputation(p,vtec_dict,re_pos,elev_rad,az_rad,gpst,rover_t,freq)
% Compute Ionospheric delay using SSR VTEC Spherical Harmonic model.

% Find the closet vtec parameters
[~, vtec_para] = closestKeyAndValue(vtec_dict, rover_t, 0, 60);
if isempty(vtec_para)
    iono_corr = NaN;
    return;
end

r_m = norm(re_pos);
[lat_ip, lon_ip, center_ang] = pierce_point_calc(deg2rad(p.lat_deg),deg2rad(p.lon_deg),az_rad,elev_rad,...
    r_m,vtec_para.height_m);

% The mean sun fixed longitude phase shifted by 2 h to the approximate 
% TEC maximum at 14:00 local time (resp. 50400 s)
epoch = mod(gpst.sow, 86400.0);
lonS = mod((lon_ip + (epoch - 50400) * pi / 43200), 2*pi);

% Compute single layer contribution
vtec = 0;
for n_deg = 0:vtec_para.num_deg
    for m_order = 0:min(n_deg, vtec_para.num_order)
         pnm = associatedLegendreFunction(n_deg, m_order, sin(lat_ip));
         a = factorial(n_deg - m_order);
         b = factorial(n_deg + m_order);
         if m_order==0
            factor = sqrt(2 * n_deg + 1);
         else
            factor = sqrt(2 * (2 * n_deg + 1) * a / b);
         end
         pnm = pnm * factor;
         Cnm_mlambda = vtec_para.cos_coeffs(n_deg+1, m_order+1) * cos(m_order * lonS);
         Snm_mlambda = vtec_para.sin_coeffs(n_deg+1, m_order+1) * sin(m_order * lonS);
         vtec = vtec + (Snm_mlambda + Cnm_mlambda) * pnm;
    end
end

if vtec < 0
    warning('Computed VTEC < 0.')
    vtec = 0;
end

mapping = 1/sin(elev_rad + center_ang);
stec = vtec * mapping;
iono_corr = stec * 40.3e16 / (freq * freq);
mapping_m = mapping * 40.3e16 / (freq * freq);
end

function [lat_f, lon_f, psi_r]= pierce_point_calc(lat_rad,lon_rad,az_rad,elev_rad,r_m,hi_m)
%%%% input:
%%% lat_0: user latitude position in radians
%%% lon_0: user longitude position in radians
%%%% azim: azimuth angle in radians
%%%% elev: elevation angle in radians
%%% r_m: user radius in meter
%%% hi_m: height in meter of ionospheric layer above the spherical Earth model

%%%% output:
%%%%lat_f: latitude of pierce point
%%%%lon_f: longitude of pierce point

% Spherical Earth.s radius of 6370 km.
Re = 6370000;

% the spherical Earth's central angle between rover position and the projection of
% the pierce point to the spherical Earth's surface.
psi_r=pi/2-elev_rad-asin(r_m*cos(elev_rad)/(Re+hi_m));
%%% eqn 1
rr1=sin(lat_rad)*cos(psi_r)+cos(lat_rad)*sin(psi_r)*cos(az_rad);

%%% latitude of the ionosphere piercing point
phi_ip=asin(rr1);

% %%%% longitude of the ionosphere piercing point
if ( (lat_rad >= 0 && tan(psi_r) * cos(az_rad) > tan(pi/2 - lat_rad))...
    || (lat_rad < 0 && -(tan(psi_r) * cos(az_rad)) > tan(pi/2 + lat_rad)) )
    lambda_ip = lon_rad + pi - asin((sin(psi_r)*sin(az_rad))/(cos(phi_ip)));
else
    lambda_ip=lon_rad+asin(sin(psi_r)*sin(az_rad)/cos(phi_ip));
end

lat_f = phi_ip;
lon_f = lambda_ip;

end

function val = associatedLegendreFunction(n, m, t)
  sum = 0.0;
  r = floor((n - m) / 2);
  for k = 0:r
    sum = sum + ((-1)^k * factorial(2*n - 2*k) / ...
                 (factorial(k) * factorial(n-k) * factorial(n-m-2*k)) ...
                 * t^(n-m-2*k));
  end
  fac = (2^-n) * ((1 - t^2)^(m/2));
  val = sum * fac;
end
