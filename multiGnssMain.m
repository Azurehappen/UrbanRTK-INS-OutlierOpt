% This code implement single frequency GNSS (GPS, GLO, GAL, BDS)
% SPS Solution using standard Iono model or Dual-freq Iono estimation
% PPP solution using real-time PPP correction.
% DGNSS performs single frequency GNSS.
% RTK Float performs single frequency GNSS with instantaneous ambiguity
% mode.
% PPP corrections:
%   SSR satellite orbit and clock correction (WHU stream)
%   Troposphere correction: IGGtrop or UNB3M
%   Ionosphere correction: SSR VTEC (CNES stream)
%   Satellite code bias (OSR): GIPP product
% clear all
% close all
%--------------------------------%
addpath('data')
addpath('parser')
addpath('timeHelper')
addpath('eph')
addpath('pos')
addpath('corr')
addpath('init')
addpath('imu')
addpath('estimation')
%--------------------------------%
% Pick the Data Number
data_num = 4;
[files, p] = dataPathLoader(data_num);
%--------------------------------%

% Read data
[p, eph, obs] = readDataFiles(p,files);
[p, obs] = loadDataAndCorr(p, files, eph, obs);
%%
p = initialParameters(p, files, eph);
p.post_mode  = p.mode_rtkfloat; % sps=Standard GNSS, ppp = PPP, dgnss = DGNSS

p.double_diff = true; % true to chose double-difference
p.elev_mark_rad  = deg2rad(10); % Set elevation cut-off
% To use Multi-GNSS and DGNSS, GPS have to be enabled.
p.enableGPS  = true; % Enable GPS
p.enableGLO  = true; % Enable GLO
p.enableGAL  = true; % Enable GAL
p.enableBDS  = true; % Enable BDS

if p.post_mode == p.mode_sps
    p.L2enable = false;        % true to enable L2 in SPS for slant iono estimation
    p.bia_type = true;         % true to disable PPP code bias
elseif p.post_mode == p.mode_ppp
    p.IGS_enable = 1;       % true to enable clock and orbit corrections from IGS
    p.tec_tmax = 15;        % Max latency of USTEC delay for use in PPP
    p.tec_tmin = 0;         % Min latency of USTEC delay for use in PPP
    p.enable_vtec = true;   % uses SSR VTEC message for Iono in PPP
end
p.state_mode = p.ins_mode; % (pos, pva, ins)_mode
if p.state_mode == p.ins_mode
    p.save_ins_pos = false;       % save INS pos to GNSS log if GNSS soln fails 
    p.save_ins_states = false;     % save INS state at INS computation times
    p.imu_data = files.imu_data;  % save IMU data to the data structure p
    p.imu_para = files.imu_para;  % save parameters to the data structure p (see readDataFiles)
elseif p.state_mode == p.pos_mode
    p.ekf_para.q_pos = 200^2;     % Pos driving noise cov for pos KF
end
p.est_mode = p.raps_ned_est;           % ekf_est, map_est, td_est, raps_ned_est

output = computeNavigationSol(p,eph,obs);

%% Plotting
ind = ~isnan(output.cost);
hor_err = output.hor_err(ind);
ver_err = abs(output.ned_err(3,ind));
meas = output.num_meas_used(ind);
compt = output.comp_time(ind);
horcov = zeros(1,length(ind));
vercov = sqrt(output.ned_cov(3,ind));
for i=1:length(ind)
    if (output.ned_cov(1,i) == 400)
        vercov(i) = NaN;
    else
        horcov(i) = norm(output.ned_cov(1:2,i));
    end
end
horcov = horcov(ind);

nonNaNCount = length(hor_err);
fprintf('Hor <= 1.0 m: %.2f%%\n', sum(hor_err <= 1.0) / nonNaNCount * 100);
fprintf('Hor <= 1.5 m: %.2f%%\n', sum(hor_err <= 1.5) / nonNaNCount * 100);
fprintf('Ver <= 3.0 m: %.2f%%\n', sum(ver_err <= 3.0) / nonNaNCount * 100);
fprintf('Hor Mean: %.2f\n', mean(hor_err));
fprintf('Hor RMS: %.2f\n', rms(hor_err));
fprintf('Hor Max: %.2f\n', max(hor_err));
fprintf('Ver Mean: %.2f\n', mean(ver_err));
fprintf('Ver RMS: %.2f\n', rms(ver_err));
fprintf('Ver Max: %.2f\n', max(ver_err));

figure
scatter(p.t,output.hor_err,'.')
xtickformat('yyyy-MM-dd HH:mm:ss')
title('Horizontal Positioning Error')
xlabel('Local time')
ylabel('Error, unit: meter');grid on

% Estimated trajectory
plotEstPosOnMap(output.pos_ecef, output.hor_err)

% Ground truth trajectory
plotEstPosOnMap(p.Grdpos.pos, zeros(1,size(p.Grdpos.pos,2)));