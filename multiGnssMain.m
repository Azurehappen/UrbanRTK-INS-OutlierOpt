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

% Initialize parameters
if exist(files.preload,'file')==2 % Check if the data already been parsed
    load(files.preload);
else
    [p, eph, obs] = readDataFiles(p,files);
    %--------------------------------%
    [p, obs] = loadDataAndCorr(p, files, eph, obs);
    save(files.preload, 'p', 'eph', 'obs');
end
%%
p = initialParameters(p, files, eph);
p.post_mode  = p.mode_rtkfloat; % sps=Standard GNSS, ppp = PPP, dgnss = DGNSS
p.imu_enable = 1;
p.double_diff = true;
p.elev_mark_rad  = deg2rad(10);
% To use Multi-GNSS and DGNSS, GPS have to be enabled.
p.enableGPS  = 1; % Enable GPS: 1 means enable, 0 means close
p.enableGLO  = 1; % Enable GLO: 1 means enable, 0 means close
p.enableGAL  = 1; % Enable GAL: 1 means enable, 0 means close
p.enableBDS  = 1; % Enable BDS: 1 means enable, 0 means close
p.L2enable = 0;
p.save_ins_pos = false;
p.est_mode = p.raps_ned_est; % map_est, td_est, raps_ned_est
p.state_mode = p.ins_mode; % POS, PVA, INS
% PPP related settings:
p.IGS_enable = 1;
p.tec_tmax = 15;
p.tec_tmin = 0;
p.enable_vtec = true;

if p.imu_enable == 1
    p.imu_data = files.imu_data;
    p.imu_para = files.imu_para;
    p.sigmaSquare_dop = 1^2;
    p.code_noise_fact_a = 300;
    p.code_noise_fact_b = 500;
end
output = computeNavigationSol(p,eph,obs);

%% Plotting
ind = ~isnan(output.cost);
hor_err = output.hor_err(ind);
ver_err = abs(output.ned_err(3,ind));
meas = output.num_meas_used(ind);
compt = output_ekf.comp_time(ind);
horcov = zeros(1,length(ind));
vercov = sqrt(output_ekf.ned_cov(3,ind));
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