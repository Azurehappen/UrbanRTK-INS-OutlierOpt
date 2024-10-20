function [files,p] = dataPathLoader(data_num)
p = struct;
switch(data_num)
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
        files.data_base = 'data/univOfTexas/asterx4_base_1hz.obs';
        % files.base_pos = [-742080.469;-5462030.972;3198339.001];
        files.base_pos = [-742080.4125;-5462031.7412;3198339.6909];
        % lever_arm = convlength([20.3504,14.44,3.6624],'in','m')';
        lever_arm = [0;0;0]; % Ground Truth used antenna position.
        [files.Grdpos.pos, files.Grdpos.vel, files.Grdpos.t,files.Grdpos.datet]...
            = readTexasTxtGt('data/univOfTexas/ground_truth.log', lever_arm);
        % Bosch IMU (Smart phone grade)
        [files.imu_data, files.imu_para] = readBoschImuData('data/univOfTexas/bosch_imu.log');
        files.imu_lever_arm = files.imu_para.R_i2b*convlength([16.1811,-6.2992,4.8425],'in','m')';
        % Lord IMU (Industy grade)
        % [files.imu_data, files.imu_para]= readLordImuData('data/univOfTexas/lord_imu.log');
        % files.imu_lever_arm = convlength([-0.461,-0.125,-0.119],'in','m')';
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
end

end 