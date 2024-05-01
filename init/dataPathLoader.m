function [files,p] = dataPathLoader(data_num)
p = struct;
switch(data_num)
    case 1
        files.eph = 'data/PPP/BRDC_20231440.rnx';
        files.obs = 'data/PPP/ppp_raps_test_1.obs';
        files.ssr = 'data/PPP/SSRA00CNE01440.23C';
        files.vtec = 'data/PPP/SSRA00CNE01440.23C';
        files.code_bias = 'data/PPP/CAS0MGXRAP_20231420_OSB.BIA';
        files.ustec_data = 'data/PPP/ustec_data/';
        files.Grdpos.pos = [-2430697.636;-4704189.056;3544329.070];
        files.Grdpos.t = NaN;
        files.preload = 'data/PPP/preload.mat';
    case 2
        files.eph = 'data/DGNSS_Moving/COM3_200712.nav';
        % files.eph = 'data/DGNSS_Moving/BRDC00WRD_194_MN.rnx';
        files.obs = 'data/DGNSS_Moving/COM3_200712.obs';
        files.ssr = [];
        files.vtec = [];
        files.code_bias = [];
        files.ustec_data = [];
        files.data_base = 'data/DGNSS_Moving/rbst0711base.obs';
        files.base_pos = [-2430697.636;-4704189.056;3544329.070];
        M  = readmatrix('data/DGNSS_Moving/groundtruth.csv');
        files.Grdpos.pos = M(:,3:5)';
        % This ground truth assumes the data is all in the same GPS week.
        % TODO: Improve to support generic time (unix_time computed from week and sow)
        files.Grdpos.t = M(:,2);
        files.preload = 'data/DGNSS_Moving/preload.mat';
    case 3
        files.eph = 'data/Tokyo/GAMG00KOR_2018353.rnx';
        % files.eph = 'data/Tokyo/brdcnav_353.18p';
        files.obs = 'data/Tokyo/rover_ublox.obs';
        files.ssr = [];
        files.vtec = [];
        files.code_bias = [];
        files.ustec_data = [];
        % files.data_base = 'data/Tokyo/CHOF00JPN_S_2018353_MO.obs';
        % files.base_pos = [-3946216.6620;3366689.8380;3698971.7310];
        files.data_base = 'data/Tokyo/base_trimble2.obs';
        files.base_pos = [-3961904.818;3348993.717;3698211.757];
        [files.Grdpos.pos, files.Grdpos.t] = readTokyoGt('data/Tokyo/reference.csv');
        files.preload = 'data/Tokyo/preload.mat';
     case 4
        % DGNSS-INS application
        p.sigmaSquare_dop = 1^2;
        p.code_noise_fact_a = 300;
        p.code_noise_fact_b = 500;
        % Septentrio receivers, use F/NAV msg (258) for GAL
        % Enable GPS, GLO, GAL, BDS
        files.eph = 'data/univOfTexas/brdm1290.19p';
        % files.eph = 'data/univOfTexas/asterx4_base.nav';
        files.obs = 'data/univOfTexas/asterx4_rover.obs';
        files.ssr = [];
        files.vtec = [];
        files.code_bias = [];
        files.ustec_data = [];
        % files.data_base = 'data/univOfTexas/SEPT1290base.obs';
        files.data_base = 'data/univOfTexas/asterx4_base_short.obs';
        files.base_pos = [-742080.469;-5462030.972;3198339.001];
        % lever_arm = [0.5169;0.3668;0.0930];
        lever_arm = [0;0;0];
        [files.Grdpos.pos, files.Grdpos.vel, files.Grdpos.t,files.Grdpos.datet]...
            = readTexasTxtGt('data/univOfTexas/ground_truth.log', lever_arm);
        % Bosch IMU
        % [files.imu_data, files.imu_para]= readImuData('data/univOfTexas/bosch_imu.log');
        % files.imu_lever_arm = [0.411,-0.160,0.123]';
        % Lord IMU
        [files.imu_data, files.imu_para]= readLordImuData('data/univOfTexas/lord_imu.log');
        files.imu_lever_arm = [-0.461,-0.125,-0.119]';
        files.preload = 'data/univOfTexas/preload.mat';
     case 5
        % PPP application
        p.sigmaSquare_dop = 2;
        p.code_noise_fact_a = 450;
        p.code_noise_fact_b = 0;
        % ublox, use I/NAV msg (517) for GAL
        % Enable GPS, GAL, BDS
        files.eph = 'data/ucr_ppp_moving/BRDC00WRD_S_2023237.rnx';
        files.obs = 'data/ucr_ppp_moving/rtk2.obs';
        files.ssr = 'data/ucr_ppp_moving/SSRA00WHU02370.23C';
        files.vtec = 'data/ucr_ppp_moving/SSRA00CNE02370.23C';
        files.code_bias = 'data/ucr_ppp_moving/CAS0MGXRAP_2023238_OSB.BIA';
        files.ustec_data = 'data/ucr_ppp_moving/ustec_data/';
        files.Grdpos = readUbxGt('data/ucr_ppp_moving/groud_truth_rtk.csv');
        files.preload = 'data/ucr_ppp_moving/preload.mat';
     case 6
        files.eph = 'data/ucr_dgnss_moving/BRDC00WRD_S_2023237.rnx';
        files.obs = 'data/ucr_dgnss_moving/rtk2.obs';
        files.ssr = [];
        files.vtec = [];
        files.code_bias = [];
        files.ustec_data = [];
        files.data_base = 'data/ucr_dgnss_moving/base_230825.23O';
        files.base_pos = [-2430697.636;-4704189.056;3544329.070];
        files.Grdpos = readUbxGt('data/ucr_dgnss_moving/groud_truth_rtk.csv');
        files.preload = 'data/ucr_dgnss_moving/preload.mat';
end

end 