function [r_ecef, v_ecef, E_k] = sat_posvel(p, eph, tm, prn, tidx,sys_type)
% Compute GNSS (GPS, GAL, BDS) satellite postion and velocity in ECEF
%
% Input:
%   p: GNSS system parameters (sysp)
%   eph: GNSS eph data
%   tm: Corrected mess. trans. time
% Output:
%   r_ecef: 3-by-1 vector of satellite position in ECEF coordinates (m)
%   v_ecef: 3-by-1 vector of satellite velocity in ECEF coordinates (m/s)
%   E_k:        eccentric anomaly
%
% Author: Wang Hu & Farzana Rahman

%-----------------%
A = eph.sqrtA(prn,tidx)^2; % semi-major axis (m)
n_0 = sqrt(p.mu/A^3); % mean motion (rad/s)
t_k = limit_tgps(posixTimeToGnssSec(tm) - eph.t_oe(prn,tidx)); % time from ephemeris reference epoch
n = n_0 + eph.Delta_n(prn,tidx); % corrected mean motion
M_k = eph.M_0(prn,tidx) + n*t_k; % mean anomaly
E_k = kepler_anomaly(eph.e(prn,tidx), M_k); % solve Kepler's equation for eccentric anomaly
sEk = sin(E_k); cEk = cos(E_k);
v_k = atan2(sqrt(1-eph.e(prn,tidx)^2)*sEk, cEk-eph.e(prn,tidx)); % true anomaly (rad)
cvk = cos(v_k); svk = sin(v_k);
Phi_k = v_k + eph.Omega(prn,tidx); % argument of latitude (rad)

% second harmonic perturbation
s2Phik = sin(2*Phi_k);
c2Phik = cos(2*Phi_k);
delta_u_k = eph.C_us(prn,tidx)*s2Phik + eph.C_uc(prn,tidx)*c2Phik; % correction to argument of latitude (rad)
delta_r_k = eph.C_rs(prn,tidx)*s2Phik + eph.C_rc(prn,tidx)*c2Phik; % correction to orbit radius in (m)
delta_i_k = eph.C_is(prn,tidx)*s2Phik + eph.C_ic(prn,tidx)*c2Phik; % correction to inclination in (rad)

% corrected argument of latitude (rad)
u_k = Phi_k + delta_u_k;
cuk = cos(u_k);
suk = sin(u_k);
c2uk = cos(2*u_k);
s2uk = sin(2*u_k);

% corrected radius (m)
r_k = A*(1 - eph.e(prn,tidx)*cEk) + delta_r_k;

% corrected inclination (rad)
i_k = eph.i_0(prn,tidx) + delta_i_k + eph.IDOT(prn,tidx)*t_k;
cik = cos(i_k);
sik = sin(i_k);

% SV positions in orbital plane after rotation thru argument of latitude (m)
x_k_prime = r_k*cuk;
y_k_prime = r_k*suk;

if strcmp('BDS',sys_type)&& (prn<=5 || prn == 18 || prn>=59)
    % corrected longitude of ascending node (rad)
    Omega_k = eph.Omega_0(prn,tidx) + eph.OmegaDot(prn,tidx)*t_k - p.OmegaDot_e*eph.t_oe(prn,tidx);
    sOmeg = sin(Omega_k);
    cOmeg = cos(Omega_k);
    % SV position in earth-fixed coordinates (ECEF)
    xg = x_k_prime*cOmeg - y_k_prime*cik*sOmeg;
    yg = x_k_prime*sOmeg + y_k_prime*cik*cOmeg;
    zg = y_k_prime*sik;
    sino = sin(p.OmegaDot_e*t_k); coso = cos(p.OmegaDot_e*t_k);
    x_k = xg*coso+yg*sino*cosd(-5)+zg*sino*sind(-5);
    y_k = -xg*sino+yg*coso*cosd(-5)+zg*coso*sind(-5);
    z_k = -yg*sind(-5)+zg*cosd(-5);
else
    % corrected longitude of ascending node (rad)
    Omega_k = eph.Omega_0(prn,tidx) + (eph.OmegaDot(prn,tidx) - p.OmegaDot_e)*t_k - p.OmegaDot_e*eph.t_oe(prn,tidx);
    sOmeg = sin(Omega_k);
    cOmeg = cos(Omega_k);
    % SV position in earth-fixed coordinates (ECEF)
    x_k = x_k_prime*cOmeg - y_k_prime*cik*sOmeg;
    y_k = x_k_prime*sOmeg + y_k_prime*cik*cOmeg;
    z_k = y_k_prime*sik;
end
r_ecef = [x_k; y_k; z_k];


% satellite velocity estimation
% Farrell (C.4.1)
E_k_dot = n/(1-eph.e(prn,tidx)*cEk);

Phi_k_dot = ((sqrt(1-eph.e(prn,tidx)^2))/(1-eph.e(prn,tidx)*cEk))* E_k_dot;

u_k_dot = (1 + 2*eph.C_us(prn,tidx)*c2Phik-2*eph.C_uc(prn,tidx)*s2Phik)*Phi_k_dot;

r_k_dot = 2*(eph.C_rs(prn,tidx)*c2Phik - eph.C_rc(prn,tidx)*s2Phik)*Phi_k_dot + A*eph.e(prn,tidx)* sEk* E_k_dot;

x_k_prime_dot = r_k_dot*cuk - r_k*suk*u_k_dot;

y_k_prime_dot = r_k_dot*suk + r_k*cuk*u_k_dot;

di_k_dt = 2*(eph.C_is(prn,tidx)*c2Phik - eph.C_ic(prn,tidx)*s2Phik)* Phi_k_dot + eph.IDOT(prn,tidx);

Omega_k_dot = eph.OmegaDot(prn,tidx) - p.OmegaDot_e;

x_k_dot = x_k_prime_dot*cOmeg - y_k_prime_dot*cik*sOmeg + y_k_prime*sik*sOmeg*di_k_dt - y_k*Omega_k_dot;

y_k_dot = x_k_prime_dot*sOmeg + y_k_prime_dot*cik*cOmeg - y_k_prime*sik*cOmeg*di_k_dt + x_k*Omega_k_dot;

z_k_dot = y_k_prime_dot*sik + y_k_prime*cik*di_k_dt;


% satellite velocity estimation
% Benjamin W. Remondi, February 2004
% http://www.ngs.noaa.gov/gps-toolbox/bc_velo/bc_velo.c

% derivative of mean anomaly (rad/s)
Mdot_k = n;

% derivative of eccentric anomaly (rad/s)
Edot_k = Mdot_k/(1.0 - eph.e(prn,tidx)*cEk);

% derivative of true anomaly (rad/s)
vdot_k = sEk*Edot_k*(1+eph.e(prn,tidx)*cvk) / (svk*(1-eph.e(prn,tidx)*cEk));

% derivative of corrected argument of latitude (rad/s)
udot_k = vdot_k + 2*(eph.C_us(prn,tidx)*c2uk - eph.C_uc(prn,tidx)*s2uk)*vdot_k;

% derivative of corrected radius (m/s)
rdot_k = A*eph.e(prn,tidx)*sEk*n/(1-eph.e(prn,tidx)*cEk) + 2*(eph.C_rs(prn,tidx)*c2uk - eph.C_rc(prn,tidx)*s2uk)*vdot_k;

% derivative of corrected inclination (rad/s)
idot_k = eph.IDOT(prn,tidx) + (eph.C_is(prn,tidx)*c2uk - eph.C_ic(prn,tidx)*s2uk)*2*vdot_k;

% derivative of positions in orbital plane
xdot_k_prime = rdot_k*cuk - y_k_prime*udot_k;
ydot_k_prime = rdot_k*sik + x_k_prime*udot_k;

% derivative of corrected longitude of ascending node (rad/s)
OmegaDot_k = eph.OmegaDot(prn,tidx) - p.OmegaDot_e;

% derivative of earth-fixed coordinates after rotation thru argument of latitude (m/s)
xdot_k = (xdot_k_prime - y_k_prime*cik*OmegaDot_k)*cOmeg ...
    - (x_k_prime*OmegaDot_k + ydot_k_prime*cik - y_k_prime*sik*idot_k)*sOmeg;
ydot_k = (xdot_k_prime - y_k_prime*cik*OmegaDot_k)*sOmeg ...
    + (x_k_prime*OmegaDot_k + ydot_k_prime*cik - y_k_prime*sik*idot_k)*cOmeg;
zdot_k = ydot_k_prime*sik + y_k_prime*cik*idot_k;

% SV velocity wrt earth in ECEF coordinates (m)
v_ecef = [xdot_k; ydot_k; zdot_k];
%     v_ecef = [x_k_dot; y_k_dot; z_k_dot];
end