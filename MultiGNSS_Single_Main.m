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
p.run_mode = 0;
p.post_mode  = p.mode_rtkfloat; % sps=Standard GNSS, ppp = PPP, dgnss = DGNSS
p.imu_enable = 1;
p.IGS_enable = 1;
p.VRS_mode = 0;
p.double_diff = true;
p.elev_mark_rad  = deg2rad(10);
% To use Multi-GNSS and DGNSS, GPS have to be enabled.
p.enableGPS  = 1; % Enable GPS: 1 means enable, 0 means close
p.enableGLO  = 1; % Enable GLO: 1 means enable, 0 means close
p.enableGAL  = 1; % Enable GAL: 1 means enable, 0 means close
p.enableBDS  = 1; % Enable BDS: 1 means enable, 0 means close
p.inval = 1; % Computation time interval
p.tec_tmax = 15;
p.tec_tmin = 0;
p.L2enable = 0;
p.enable_vtec = true;
p.save_ins_pos = false;
p.est_mode = p.raps_ned_est; % ekf_est, map_est, td_est, raps_ned_est
p.state_mode = p.ins_mode; % POS, PVA, INS
% obs = p.obs_b;
% p.Grdpos.pos = [-742080.469;-5462030.972;3198339.001];
% p.Grdpos.t = NaN;
if p.state_mode == p.pos_mode
    p.ekf_para.q_pos = 200^2;
end
if p.imu_enable == 1
    p.imu_data = files.imu_data;
    p.imu_para = files.imu_para;
    p.sigmaSquare_dop = 1^2;
    p.code_noise_fact_a = 300;
    p.code_noise_fact_b = 500;
end
output = compute_gnss_ecef(p,eph,obs);

%%
nonNaNCount = sum(~isnan(output.hor_err));
percentage = sum(output.hor_err < 1.0) / nonNaNCount * 100
percentage = sum(output.hor_err < 1.5) / nonNaNCount * 100
percentage = sum(output.err < 3.0) / nonNaNCount * 100

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
scatter(p.t,output.open_sky,'.')
xtickformat('yyyy-MM-dd HH:mm:ss')
title('Open sky indicator')
xlabel('Local time')
ylabel('Indicator');grid on

total = output.sv_num_GPS + output.sv_num_GLO + output.sv_num_GAL + output.sv_num_BDS;
figure
scatter(p.t,total,'.')
title('total satellites been used')
% legend('GPS','GAL','BDS','Total')
xlabel('Receiver time using GPS second')
ylabel('Distance, unit: meter');grid on
% %
figure
scatter(p.t,output.rover_clk/p.c,'.')
title('Local bias')
xlabel('Receiver time using GPS second');
ylabel('Clock bias, seconds');grid on

if p.enableGLO
    figure
    scatter(p.t,output.clk_glo-output.clk_gps,'.')
    title('ISB GLO')
    xlabel('Receiver time using GPS second');
    ylabel('GPS-GLO ISB, meter');grid on
end

if p.enableGAL == 1
    figure
    scatter(p.t,output.clk_gal-output.clk_gps,'.')
    title('ISB GAL')
    xlabel('Receiver time using GPS second');
    ylabel('GPS-GAL ISB, meter');grid on
end

if p.enableBDS == 1
    figure
    scatter(p.t,output.clk_bds-output.clk_gps,'.')
    title('ISB BDS')
    xlabel('Receiver time using GPS second');
    ylabel('GPS-BDS ISB, meter');grid on
end

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

if p.est_mode == p.raps_est || p.est_mode == p.raps_ned_est
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

if p.state_mode == p.pva_mode
    figure
    subplot(311)
    grid on
    plot(p.t,output.vel_ned(1,:),'.')
    subplot(312)
    plot(p.t,output.vel_ned(2,:),'.')
    subplot(313)
    plot(p.t,output.vel_ned(3,:),'.')
    title('Plot for Velocity')
    xlabel('Local time')
    ylabel('NED Velocity, m/s')
    if isfield(p.Grdpos, 'vel')
        figure
        subplot(311)
        plot(p.t,output.vel_ned_err(1,:),'.')
        subplot(312)
        plot(p.t,output.vel_ned_err(2,:),'.')
        subplot(313)
        plot(p.t,output.vel_ned_err(3,:),'.')
        grid on
        legend('North', 'East', 'Vertical');
        title('Error Plot for Velocity')
        xlabel('Local time')
        ylabel('NED Velocity Error, m/s')
        figure
        subplot(311)
        plot(p.t,output.vel_ned(1,:),'.')
        hold on
        plot(p.Grdpos.datet,p.Grdpos.vel(1,:))
        legend('North Exp', 'North Gt');
        title('Rover speed in NED')
        subplot(312)
        plot(p.t,output.vel_ned(2,:),'.')
        hold on
        plot(p.Grdpos.datet,p.Grdpos.vel(2,:))
        legend('East Exp', 'East Gt');
        subplot(313)
        plot(p.t,-output.vel_ned(3,:),'.')
        hold on
        plot(p.Grdpos.datet,p.Grdpos.vel(3,:))
        legend('Vertivcal Exp', 'Vertical Gt');
    end
end

plotEstPosOnMap(output.pos_ecef, output.hor_err)

plotEstPosOnMap(p.Grdpos.pos, zeros(1,size(p.Grdpos.pos,2)));
% A = zeros(1, length(p.Grdpos.pos));
% plotEstPosOnMap(p.Grdpos.pos, A)
% figure
% subplot(311)
% scatter(p.t,output.ned_err(1,:),'.')
% title('North Error in NED');grid on;
% subplot(312)
% scatter(p.t,output.ned_err(2,:),'.')
% title('East Error in NED');grid on;
% subplot(313)
% scatter(p.t,output.ned_err(3,:),'.')
% title('Down Error in NED')
% xlabel('Receiver time using GPS second');grid on;
% figure
% for i=1:32
%     scatter(output.gpst,output.res_GPS(i,:),'.')
%     hold on
% end
% hold off
% grid on
% title('GPS residual')
% xlabel('Receiver time using GPS second');
% ylabel('Residual, unit: meter');
% figure
% scatter(output.gpst,output.sv_num_GPS,'.')
% title('Number of GPS satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on
%
% figure
% scatter(output.gpst,output.sv_num_GLO,'.')
% title('Number of GLO satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on
%
% figure
% scatter(output.gpst,output.sv_num_GAL,'.')
% title('Number of GAL satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on
%
% figure
% scatter(output.gpst,output.sv_num_BDS,'.')
% title('Number of BDS satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on