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
addpath('time_compute')
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
p.est_mode = p.raps_ned_est; % ekf_est, map_est, td_est, raps_ned_est
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
output = compute_gnss_ecef(p,eph,obs);

%%
figure
scatter(p.t,output.err,'.')
xtickformat('yyyy-MM-dd HH:mm:ss')
title('ECEF positioning error')
xlabel('Local time')
ylabel('Error, unit: meter');grid on

figure
scatter(p.t,output.hor_err,'.')
xtickformat('yyyy-MM-dd HH:mm:ss')
title('Horizontal positioning error')
xlabel('Local time')
ylabel('Error, unit: meter');grid on

figure
subplot(311)
scatter(p.t,sqrt(output.ned_cov(1,:)),'.')
title('NED Pos Covariance')
ylabel('Cov N, meter');grid on
subplot(312)
scatter(p.t,sqrt(output.ned_cov(2,:)),'.')
ylabel('Cov E, meter');grid on
subplot(313)
scatter(p.t,sqrt(output.ned_cov(3,:)),'.')
ylabel('Cov D, meter');grid on
title('NED Pos Covariance')
xlabel('Receiver time using GPS second');

figure
subplot(311)
scatter(p.t,sqrt(output.state_cov(1,:)),'.')
ylabel('Cov x, meter');grid on
subplot(312)
scatter(p.t,sqrt(output.state_cov(2,:)),'.')
ylabel('Cov y, meter');grid on
subplot(313)
scatter(p.t,sqrt(output.state_cov(3,:)),'.')
ylabel('Cov z, meter');grid on
title('ECEF Pos Covariance')
xlabel('Receiver time using GPS second');

if p.est_mode == p.raps_ned_est
    figure
    scatter(p.t,total - output.raps_num_sat,'.')
    title('No. of satellites been removed')
    xlabel('Receiver time using GPS second')
    ylabel('Distance, unit: meter');grid on

    figure
    subplot(311)
    scatter(p.t,sqrt(output.pos_info_ned(1,:)),'.')
    hold on
    yline(sqrt(1/p.raps.poshor_cov_spec))
    ylabel('Infor N, meter');grid on
    subplot(312)
    scatter(p.t,sqrt(output.pos_info_ned(2,:)),'.')
    hold on
    yline(sqrt(1/p.raps.poshor_cov_spec))
    ylabel('Infor E, meter');grid on
    subplot(313)
    scatter(p.t,sqrt(output.pos_info_ned(3,:)),'.')
    hold on
    yline(sqrt(1/p.raps.posver_cov_spec))
    ylabel('Infor D, meter');grid on
    title('NED Pos Information')
    xlabel('Receiver time using GPS second');
end

% Estimated trajectory
plotEstPosOnMap(output.pos_ecef, output.hor_err)

% Ground truth trajectory
plotEstPosOnMap(p.Grdpos.pos, zeros(1,size(p.Grdpos.pos,2)));