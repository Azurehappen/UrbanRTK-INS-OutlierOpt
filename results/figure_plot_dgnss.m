load('EKF-INS-RTK.mat')
output_ekf = output;
p_ekf = p;

load('TD-INS-RTK.mat')
output_td = output;
p_td = p;

load('RAPS-INS-RTK.mat')
output_raps = output;
p_raps = p;

load('RAPS-PVA-RTK.mat')
output_rapspva = output;
p_rapspva = p;
%%
ind = ~isnan(output_ekf.cost);
ekf_hor_err = output_ekf.hor_err(ind);
ekf_ver_err = abs(output_ekf.ned_err(3,ind));
ekf_meas = output_ekf.num_meas_used(ind);
ekf_compt = output_ekf.comp_time(ind);
ekf_gdop = output_ekf.GDOP(ind);
ekf_horcov = zeros(1,length(ind));
ekf_vercov = sqrt(output_ekf.ned_cov(3,ind));
for i=1:length(ind)
    if (output_ekf.ned_cov(1,i) == 400)
        ekf_vercov(i) = NaN;
    else
        ekf_horcov(i) = norm(output_ekf.ned_cov(1:2,i));
    end
end
ekf_horcov = ekf_horcov(ind);

ind = ~isnan(output_td.cost);
td_hor_err = output_td.hor_err(ind);
td_ver_err = abs(output_td.ned_err(3,ind));
td_meas = output_td.num_meas_used(ind);
td_compt = output_td.comp_time(ind);
td_gdop = output_td.GDOP(ind);
td_horcov = zeros(1,length(ind));
td_vercov = sqrt(output_td.ned_cov(3,ind));
for i=1:length(ind)
    if (output_ekf.ned_cov(1,i) == 400)
        td_vercov(i) = NaN;
    else
        td_horcov(i) = norm(output_td.ned_cov(1:2,i));
    end
end
td_horcov = td_horcov(ind);

rapspva_hor_err = output_rapspva.hor_err(ind);
rapspva_ver_err = abs(output_rapspva.ned_err(3,ind));
rapspva_risk = output_rapspva.pos_risk(ind);
rapspva_meas = output_rapspva.num_meas_used(ind);
rapspva_penalty = output_rapspva.raps_penalty(ind);
rapspva_compt = output_rapspva.comp_time(ind);
rapspva_gdop = output_rapspva.GDOP(ind);
rapspva_horcov = zeros(1,length(ind));
rapspva_vercov = sqrt(output_rapspva.ned_cov(3,ind));
for i=1:length(ind)
    if (output_rapspva.ned_cov(1,i) == 400)
        rapspva_vercov(i) = NaN;
    else
        rapspva_horcov(i) = norm(output_rapspva.ned_cov(1:2,i));
    end
end
rapspva_horcov = rapspva_horcov(ind);

raps_hor_err = output_raps.hor_err(ind);
raps_ver_err = abs(output_raps.ned_err(3,ind));
raps_risk = output_raps.pos_risk(ind);
raps_meas = output_raps.num_meas_used(ind);
raps_penalty = output_raps.raps_penalty(ind);
raps_compt = output_raps.comp_time(ind);
raps_gdop = output_raps.GDOP(ind);
raps_horcov = zeros(1,length(ind));
raps_vercov = sqrt(output_raps.ned_cov(3,ind));
for i=1:length(ind)
    if (output_raps.ned_cov(1,i) == 400)
        raps_vercov(i) = NaN;
    else
        raps_horcov(i) = norm(output_raps.ned_cov(1:2,i));
    end
end
raps_horcov = raps_horcov(ind);

nonNaNCount = length(ekf_hor_err);
fprintf('EKF RTK Hor <= 1.0 m: %.2f%%\n', sum(ekf_hor_err <= 1.0) / nonNaNCount * 100);
fprintf('EKF RTK Hor <= 1.5 m: %.2f%%\n', sum(ekf_hor_err <= 1.5) / nonNaNCount * 100);
fprintf('EKF RTK Ver <= 3.0 m: %.2f%%\n', sum(ekf_ver_err <= 3.0) / nonNaNCount * 100);
fprintf('EKF RTK Hor Mean: %.2f\n', mean(ekf_hor_err));
fprintf('EKF RTK Hor RMS: %.2f\n', rms(ekf_hor_err));
fprintf('EKF RTK Hor Max: %.2f\n', max(ekf_hor_err));
fprintf('EKF RTK Ver Mean: %.2f\n', mean(ekf_ver_err));
fprintf('EKF RTK Ver RMS: %.2f\n', rms(ekf_ver_err));
fprintf('EKF RTK Ver Max: %.2f\n', max(ekf_ver_err));

nonNaNCount = length(td_hor_err);
fprintf('TD RTK Hor <= 1.0 m: %.2f%%\n', sum(td_hor_err <= 1.0) / nonNaNCount * 100);
fprintf('TD RTK Hor <= 1.5 m: %.2f%%\n', sum(td_hor_err <= 1.5) / nonNaNCount * 100);
fprintf('TD RTK Ver <= 3.0 m: %.2f%%\n', sum(td_ver_err <= 3.0) / nonNaNCount * 100);
fprintf('TD RTK Hor Mean: %.2f\n', mean(td_hor_err));
fprintf('TD RTK Hor RMS: %.2f\n', rms(td_hor_err));
fprintf('TD RTK Hor Max: %.2f\n', max(td_hor_err));
fprintf('TD RTK Ver Mean: %.2f\n', mean(td_ver_err));
fprintf('TD RTK Ver RMS: %.2f\n', rms(td_ver_err));
fprintf('TD RTK Ver Max: %.2f\n', max(td_ver_err));

nonNaNCount = length(raps_hor_err);
fprintf('RAPS RTK Hor <= 1.0 m: %.2f%%\n', sum(raps_hor_err <= 1.0) / nonNaNCount * 100);
fprintf('RAPS RTK Hor <= 1.5 m: %.2f%%\n', sum(raps_hor_err <= 1.5) / nonNaNCount * 100);
fprintf('RAPS RTK Ver <= 3.0 m: %.2f%%\n', sum(raps_ver_err <= 3.0) / nonNaNCount * 100);
fprintf('RAPS RTK Hor Mean: %.2f\n', mean(raps_hor_err));
fprintf('RAPS RTK Hor RMS: %.2f\n', rms(raps_hor_err));
fprintf('RAPS RTK Hor Max: %.2f\n', max(raps_hor_err));
fprintf('RAPS RTK Ver Mean: %.2f\n', mean(raps_ver_err));
fprintf('RAPS RTK Ver RMS: %.2f\n', rms(raps_ver_err));
fprintf('RAPS RTK Ver Max: %.2f\n', max(raps_ver_err));

nonNaNCount = length(rapspva_hor_err);
fprintf('RAPS PVA RTK Hor <= 1.0 m: %.2f%%\n', sum(rapspva_hor_err <= 1.0) / nonNaNCount * 100);
fprintf('RAPS PVA RTK Hor <= 1.5 m: %.2f%%\n', sum(rapspva_hor_err <= 1.5) / nonNaNCount * 100);
fprintf('RAPS PVA RTK Ver <= 3.0 m: %.2f%%\n', sum(rapspva_ver_err <= 3.0) / nonNaNCount * 100);
fprintf('RAPS PVA RTK Hor Mean: %.2f\n', mean(rapspva_hor_err));
fprintf('RAPS PVA RTK Hor RMS: %.2f\n', rms(rapspva_hor_err));
fprintf('RAPS PVA RTK Hor Max: %.2f\n', max(rapspva_hor_err));
fprintf('RAPS PVA RTK Ver Mean: %.2f\n', mean(rapspva_ver_err));
fprintf('RAPS PVA RTK Ver RMS: %.2f\n', rms(rapspva_ver_err));
fprintf('RAPS PVA RTK Ver Max: %.2f\n', max(rapspva_ver_err));

plotEstPosOnMap(output_td.pos_ecef, output_td.hor_err)
plotEstPosOnMap(output_raps.pos_ecef, output_raps.hor_err)

%%
purple = [0.4940, 0.1840, 0.5560];
blue = [0, 0.4470, 0.7410];
red = [0.6350, 0.0780, 0.1840];
green = [0.4660, 0.6740, 0.1880];
orange = [0.9290, 0.6940, 0.1250];