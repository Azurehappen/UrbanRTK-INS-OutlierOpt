function [p, eph, obs] = initialParameters(p, files, eph)

p.mode_sps = 0;
p.mode_ppp = 1;
p.mode_dgnss = 2;
p.mode_rtkfloat = 3;
p.mode_rtkfix = 4;

%----------------------%
% load('data/DCB_GLO.mat');
% p.icb_glo = DCB_P1C1;
% Setting
p.inval = 1; % Computation time interval
p.run_mode=1; %%%%% 0=real-time; 1=post processing
p.post_mode = 1;%%%% 0=Standard GNSS, 1 = PPP, 2= DGNSS
p.iono_map=0; %%%%% 0=USTEC; 1=IGS ionex
p.ins=1;    %%%%%%% 0=no imu data; 1=use simulated/real imu data
p.imu_freq=200; %%%%%% frequency rate for simulated/real IMU data
if ~isempty(p.obs_b)
    if ~isfield(files, 'base_pos')
        error('Base obs available, but base pos is not provided.')
    end
    p.P_base = files.base_pos; % Base station position
end
p.Grdpos = files.Grdpos; % Ground truth of rover
p.freq = 1; %%%% 1 = single frequency, 2 = dual frequency
p.Ek0 = 0; % Initial condition of Ek
p.enableGPS = 1; % Enable GPS: 1 means enable, 0 means close
p.enableGLO = 0; % Enable GLO: 1 means enable, 0 means close
p.enableGAL = 1; % Enable GAL: 1 means enable, 0 means close
p.enableBDS = 1; % Enable BDS: 1 means enable, 0 means close
p.IGS_enable = 1; % Enable IGS correction: 1 means enable, 0 means close

p.L1freq = 1575.42e6; % L1 frequency (Hz)
p.L2freq = 1227.6e6; % L2 frequency (Hz)
p.L1glo = 1602.0e6; % GLO L1 frequency (Hz)
p.L2glo = 1246.0e6; % GLO L2 frequency (Hz)
p.E1freq = 1575.42e6; % E1 frequency (Hz)
p.E6freq = 1278.75e6; % E6 frequency (Hz)
p.E5freq = 1191.795e6; % E5 frequency (Hz)
p.E5afreq = 1176.45e6; % E5a frequency (Hz)
p.E5bfreq = 1207.14e6; % E5b frequency (Hz)
p.B1freq = 1561.098e6; % B1I frequency (Hz)
p.B2afreq = 1207.14e6; % B2b frequency (currently not supported)
p.state0 = [0;0;0;0]; % Initial state vector at first iteration, maybe changed in other functions
%                             [x;y;z;clk_bias]
p.mk=0;
p.lat = 0; p.lon = 0; p.h_r = 0; % Initialize latitude, longitude and height
p.eph_base = [];
p.obs_base = [];
% Common parameters
p.c = 2.99792458e8; % Speed of light (m/s)
p.pi = pi; % Archimedes' constant
if eph.LeapSeconds~=0
    p.gpstime_ahead = eph.LeapSeconds; % The difference between GPS time and UTC, GPS has 18 seconds ahead of UTC
else
    % In some case, there is no leap seconds in eph data.
    % This constant need to be manually changed if Leap seconds change.
    p.gpstime_ahead = 18;
end
p.Re        = 6378136.3; %%%% radius of earth in m  
p.h_iono    = 350000;    %%%% height of ionosphere of maximum TEC
% Earth model
p.g = 9.80665;                    % standard gravity
p.omge = 7.2921151467e-5;     % earth angular velocity (rad/s)
p.omg_iee_vec = [0;0;p.omge]; % earth angular vector
p.GM = 6.6743e-11*5.972e24;
p.Omega_iee_mat = vectorSkewSymMat(p.omg_iee_vec);
p.a = 6378137.0;                  % semi-major axis length, meters
p.b = 6356752.31424518;           % semi-minor axis length, meters
p.f = 1.0/298.257223563;          % flatness of ellipsoid
p.e = sqrt((p.a^2-p.b^2)/(p.a^2));      % first eccentricity of ellipsoid
p.ep = sqrt((p.a^2-p.b^2)/(p.b^2));     % second eccentricity of ellipsoid
%-------------------------------------------------------------------------%
% Treshold
p.elev_mark_rad = deg2rad(10); % Elevation treshold
p.open_sky_elev_rad = deg2rad(5); % Open sky check minimum elevation requirement
p.sig_strg = 20; % Signal strength treshold
p.satdelta = 1e-6; % In function 'eph2pos', the threshold for sat pos convergence
p.min_sv = 5; % the minimum of the number of satellites
p.LSthrsh = 1e-8; % threshold of delta_x in LS solver
% Iteration Number
p.NsatEk = 20; % In function 'eph2pos', the maxinum of iteration for Ek convergence
p.Nls = 20; % Least square
p.GDOP_mark = 30;
p.pivot_elev = deg2rad(45);
%-------------------------------------------------------------------------%
% Time Synchronization
% GPS starting at 1980-1-6 00:00:00
% Galileo starting at 1999-8-22 00:00:13
% BDS starting at 2006-1-1 00:00:00
% Leap seconds between 1980 and 1999 is 13s, Hence GPS seconds = GAL seconds
p.gal.lps_gps = 0;
% Leap seconds between 1980 and 2060 is 14s, Hence GPS seconds = BDS seconds + 14
p.bds.lps_gps = 14;
% GPS has 18 seconds ahead of UTC, GLONASS is synchronized to UCT
% Hence GPS seconds = GLO seconds + p.gpstime_ahead
p.glo.lps_gps = p.gpstime_ahead;

%---------------------------------------%  
% GPS Constants
% Reference 'https://www.gps.gov/technical/icwg/IS-GPS-200H.pdf'
%---------------------------------------%
p.gps.sys_num = 1;
p.gps.mu = 3.986005e+14; % WGS84 value of the earth's gravitational constant for GPS user (m^3/s^2)
p.gps.OmegaDot_e = 7.2921151467e-5; % Earth rotation rate (rad/s)
p.gps.F = -4.442807633e-10; % IS-GPS-200H page 96
p.gps.message_duration = 7200; % Maximum time difference between obs data and eph message
p.gps.l1_wavelength = p.c/p.L1freq;
p.gps.l2_wavelength = p.c/p.L2freq;
% GLONASS Constants
% Reference ''
%---------------------------------------%
p.glo.sys_num = 2;
p.glo.df1 = 0.56250e6; % GLONASS G1 bias frequency (Hz/n)
p.glo.mu = 3.986004418e+14; % Geocentric gravitational constant (m^3/s^2)
p.glo.OmegaDot_e = 7.2921151467e-5; % Earth's rotation rate
p.glo.F = -2*sqrt(p.glo.mu)/(p.c^2); % 
p.glo.a_e = 6378136; %  semi-major (equatorial) axis of the PZ-90 Earth’s ellipsoid
p.glo.C_20 = 1082625.75e-9; %  second degree zonal coefficient of normal potential
p.glo.message_duration = 1800; %54000;
% Galileo Constants
% Reference 'https://www.gsc-europa.eu/sites/default/files/sites/all/files/Galileo-OS-SIS-ICD.pdf'
%---------------------------------------%
p.gal.sys_num = 3;
p.gal.mu = 3.986004418e+14; % Geocentric gravitational constant (m^3/s^2)
p.gal.OmegaDot_e = 7.2921151467e-5; % mean angular velocity of the Earth
p.gal.F = -4.442807309e-10; % Galileo-OS-SIS-ICD page 58
p.gal.message_duration = 7200;
p.gal.e1_wavelength = p.c/p.E1freq;
% BeiDou Constants
% Reference 'http://en.beidou.gov.cn/SYSTEMS/ICD/201902/P020190227702348791891.pdf'
%---------------------------------------%
p.bds.sys_num = 4;
p.bds.mu = 3.986004418e+14; % Geocentric gravitational constant (m^3/s^2)
p.bds.OmegaDot_e = 7.2921150e-5; % Earth's rotation rate
p.bds.F = -2*sqrt(p.bds.mu)/(p.c^2); % BeiDou-ICD page 58
p.bds.message_duration = 3600;
p.bds.e1_wavelength = p.c/p.B1freq;
% Measurement selection
p.select = 0;
p.bia_type = 0; % 0 means using CNE code bias
%---------------------------------------%
p.GPS_C1C = 1;p.GPS_C1W = 2;p.GPS_C2L = 3;p.GPS_C2W = 4;
p.GLO_C1C = 1;p.GLO_C1P = 2;p.GLO_C2C = 3;p.GLO_C2P = 4;
p.GAL_C1X = 1;p.GAL_C7X = 2;
p.BDS_C2I = 1;p.BDS_C7I = 2;

% ISB model:
% continuous time: dx(t) = u x(t) + w (See 4.6.5 in Farrell's book)
% dP(t) = 2 u P(t) + Q
% ISB is a constant, therefore, dP(t) = 0.
% Discrete time:
% x(k+1) = Phi * x(k)
% P(k+1) = Phi * P(k) * Phi' + Q
p.ekf_para.isb_cov = 1;
p.ekf_para.u_isb = -1e-4;
p.ekf_para.q_isb = -2 * p.ekf_para.u_isb * p.ekf_para.isb_cov;
%p.ekf_para.q_isb = 5^2;

%---------------------------------------%
% PPP parameter
p.ppp.sigma_cme = 0.1; % PPP correction RMS except Iono (m)
p.ppp.sigma_iono = 6.15; % Iono correction RMS (TECU)

%---------------------------------------%
% estimation mode
p.ekf_est = 1;
p.map_est = 2;
p.td_est = 3;
p.raps_est = 4;
p.raps_ned_est = 5;
p.est_mode = p.ekf_est;

%---------------------------------------%
% EKF state mode
p.pos_mode = 1;
p.pva_mode = 2;
p.ins_mode = 3;
p.state_mode = p.pva_mode;
p.modeToNumUserErrStates = dictionary;
p.modeToNumUserErrStates(p.pos_mode) = 3;
p.modeToNumUserErrStates(p.pva_mode) = 9;
p.modeToNumUserErrStates(p.ins_mode) = 6+3+6;

%---------------------------------------%
% EKF parameters
p.ekf_para.q_pos = 100^2; % Only be applied to Pos model.
p.ekf_para.q_vel = 50^2;
p.ekf_para.q_accHor = 4^2;
p.ekf_para.q_accVer = 1^2;
p.ekf_para.q_clkDrift = 5.0^2;

p.td_lambda = 2;

%---------------------------------------%
% RAPS parameters
p.raps.clk_cov_spec = (500^2)*0.05;
p.raps.dclk_cov_spec = (50^2)*0.05;
p.raps.isb_cov_spec = (10^2)*0.05;
p.raps.va_cov_spec = (10^2)*0.05; 
% Horizontal: alpha = 1.5, beta = 0.32; 
% Vertical: alpha = 5, beta = 0.4.
p.raps.poshor_cov_spec = 0.72; % 1.5^2*0.32
p.raps.posver_cov_spec = 2.88; % 3^2*0.32
p.raps.velhor_cov_spec = 0.72/2;
p.raps.velver_cov_spec = 2.88/2;
p.raps_penalty = NaN;
p.check_prior = false;

p.num_sats_window = NaN(1,10);
p.GDOP = NaN;

if isfield(files, 'imu_lever_arm')
    p.imu_lever_arm = files.imu_lever_arm;
end

end
