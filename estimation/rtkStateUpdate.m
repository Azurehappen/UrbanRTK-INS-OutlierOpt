function [p,estState,res] = rtkStateUpdate(p, cpt, dt)
% RTK Instantaneous
%-------------------%
% Initialize
estState.clock_sys = dictionary;
estState.clock_sys(p.gps.sys_num) = NaN;
estState.clock_sys(p.glo.sys_num) = NaN;
estState.clock_sys(p.gal.sys_num) = NaN;
estState.clock_sys(p.bds.sys_num) = NaN;
x_minus = p.state0;
num_user_states = p.modeToNumUserStates(p.state_mode);
% [H_clk,x_clk] = formClkStatesAndH(cpt.num_sv);
% if length(x_clk) + num_user_states + 1 ~= length(p.state0)
%     error('current No. of sys does not match the previous epoch');
% end
x_clk = zeros(p.enableGPS+p.enableGLO+p.enableGAL+p.enableBDS,1);
H_clk = zeros(sum(cpt.num_sv),length(x_clk));
row_ind = 0;
col_ind = 1;
for i = 1:length(cpt.num_sv)
    if cpt.num_sv(i) == 0
        if i ~= 2
            col_ind = col_ind + 1;
        end
        continue;
    end
    H_clk(row_ind+1:row_ind+cpt.num_sv(i),col_ind)=1;
    row_ind = row_ind + cpt.num_sv(i);
    col_ind = col_ind + 1;
end
%------------------%
y_rho = cpt.corr_range;
p.num_sats_window = [p.num_sats_window(2:length(p.num_sats_window)), length(y_rho)];
y_dop = [];
num_rho = length(y_rho); % The number of measurement
num_dop = 0;
if p.state_mode == p.pva_mode
    y_dop = cpt.doppler;
    num_dop = length(y_dop);
end

y_phi = [];
num_phi = 0;
T_dd = [];
num_sv_rtk = cpt.num_sv;
if p.post_mode == p.mode_rtkfloat
    y_phi = cpt.phase_m;
    num_phi = length(y_phi);
    % Find the pivot satellite
    sysnumToPivotPrn = dictionary;
    T_dd = eye(num_phi); % Phase double-diff transformation matrix
    to_delete = [];
    st_ind = 1;
    lambda = cpt.wavelength;
    for i = 1:length(cpt.num_sv)
        ind_i = find(cpt.svprn_mark == i);
        if isempty(ind_i)
            continue;
        end
        elev_i = cpt.elev(ind_i);
        prn_i = cpt.prn_record(ind_i);
        [elev_max, ind_pivot] = max(elev_i);
        if elev_max < p.pivot_elev
            % No suitable pivot sat, remove this system
            to_delete = [to_delete,st_ind:st_ind+cpt.num_sv(i)-1];
            T_dd(:,st_ind:st_ind+cpt.num_sv(i)-1) = 0;
            num_sv_rtk(i) = 0;
            continue;
        else
            to_delete = [to_delete,st_ind+ind_pivot(1)-1];
            T_dd(st_ind:st_ind+cpt.num_sv(i)-1,st_ind+ind_pivot(1)-1) = -1;
        end
        st_ind = st_ind + cpt.num_sv(i);
        sysnumToPivotPrn(i) = prn_i(ind_pivot(1));
    end
    T_dd(to_delete,:) = [];
    lambda(to_delete) = [];
    x_minus = [x_minus;zeros(length(lambda),1)];
    H_phi = zeros(num_phi, length(x_minus));
    if sum(num_sv_rtk) < 6
        p.i
        warning('Unconfident number of sat to perform RTK');
        estState.pos = [];
        res = NaN;
        return;
    end
end

H_rho = zeros(num_rho, length(x_minus));
H_rho(:, num_user_states+1:num_user_states+length(x_clk)) = H_clk;
H_dop = [];
if p.state_mode == p.pva_mode
    H_dop = zeros(num_dop, length(x_minus));
    H_dop(:, num_user_states+length(x_clk)+1) = 1;
end

s_pos_ecef = cpt.s_pos_ecef;
if p.post_mode == p.mode_ppp && p.IGS_enable == 1
    s_pos_ecef = cpt.s_pos_prc;
end
s_v_ecef = cpt.sat_v_Rcorr;

Range = zeros(num_rho,1);
geo_range_rho = zeros(num_rho,1);
geo_range_phi = zeros(num_phi,1);
dop_sat = zeros(num_dop,1);
dd_phi = T_dd * y_phi; % double-diff phase
meas = [y_rho;T_dd * y_phi;y_dop];
R_dop = eye(length(y_dop));

[R_rho, R_phi] = constructMeasNoise(p, cpt, 1);
dd_rho = T_dd * y_rho; % double-diff pseudorange
Rt_phi = T_dd * R_phi * T_dd';
% dd_amb = (dd_phi - dd_rho)./lambda;
dd_amb = (dd_phi - dd_rho);
% cov_dd_amb = T_dd * diag((diag(R_phi)+diag(R_rho))./(cpt.wavelength.^2)) * T_dd';
cov_dd_amb = T_dd * diag((diag(R_phi)+diag(R_rho))) * T_dd'+1000^2;
x_minus(length(x_minus)-length(lambda)+1:end) = dd_amb;
cov_minus = blkdiag(p.state_cov, cov_dd_amb);
R = blkdiag(R_rho,Rt_phi,R_dop);

for j=1:num_rho
    Range(j)=norm(s_pos_ecef(:,j)-x_minus(1:3));
    los = (x_minus(1:3)-s_pos_ecef(:,j))'/Range(j)+...
        [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0];
    H_rho(j,1:3)=los;
    geo_range_rho(j) = Range(j)+sagnac(p,s_pos_ecef(:,j),x_minus(1:3));
    if p.state_mode == p.pva_mode
        H_dop(j,4:6)=los;
        range_r = norm(cpt.sat_pos_Rcorr(:,j)-x_minus(1:3));
        los_r = (x_minus(1:3)-cpt.sat_pos_Rcorr(:,j))'/range_r;
        dop_sat(j) = -los_r*s_v_ecef(:,j);
    end
    if p.post_mode == p.mode_rtkfloat
        H_phi(j,1:3)=los;
        geo_range_phi(j) = geo_range_rho(j);
    end
end
if p.post_mode == p.mode_rtkfloat
    geo_range_phi = T_dd * geo_range_phi;
    H_phi_dd = T_dd * H_phi;
    % H_phi_dd(:, length(x_minus)-length(lambda)+1:end) = diag(lambda);
    H_phi_dd(:, length(x_minus)-length(lambda)+1:end) = diag(ones(length(lambda),1));
    H_os = [H_rho;H_phi_dd;H_dop];
    geo_range = [geo_range_rho;geo_range_phi;dop_sat];
else
    H_os = [H_rho;H_dop];
    geo_range = [geo_range_rho;dop_sat];
end

% measurement residual
res_all = meas - geo_range - H_os(:,4:end)*x_minus(4:end);
[p.state0, p.state_cov, flag] = checkClockReset(p, p.state0, p.state_cov, ...
    num_user_states, res_all(1:num_rho), cpt);
if flag == true
    res_all = meas - geo_range - H_os(:,4:end)*x_minus(4:end);
    cov_minus = blkdiag(p.state_cov, cov_dd_amb);
    x_minus(1:length(x_minus)-length(lambda)) = p.state0;
end

% y - f(x0) = H (x - x0);
% zk = res_all + H_os * x_minus;
switch p.est_mode
    case p.ekf_est
        [x_plus, cov_plus] = ekfUpdate(x_minus, p.state_cov, res_all, H_os, R);
    % case p.td_est
    %     [x_plus, cov_plus] = ekfUpdate(x_minus, p.state_cov, res_all, H_os, R);
    case p.map_est
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla');
        R_pva = [R_e2g, zeros(3,6);
            zeros(3,3), R_e2g, zeros(3,3);
            zeros(3,6), R_e2g];
        Rot_e2g = [R_pva, zeros(9,length(x_minus)-9);
            zeros(length(x_minus)-9, 9), eye(length(x_minus)-9)];
        H_os = H_os * Rot_e2g';
        cov_prior = Rot_e2g' * p.state_cov * Rot_e2g;
        [x_plus,cov_plus,p.infor_ned,p.augcost] = ...
            mapUpdate(ones(num,1),x_minus,cov_prior,res_all,H_os,R,Rot_e2g);
        p.num_meas_used = num;
    case p.td_est
        % [flag_rapid,p.num_sats_window] = checkRapidNumSatChange(p.num_sats_window, sum(cpt.num_sv~=0));
        % if p.state_mode == p.pva_mode && flag_rapid == true
        %     p.state_cov = zeros(size(p.state_cov));
        %     p.state_cov(1:end-1,1:end-1) = diag(150^2*ones(1,length(x_minus)-1));
        %     p.state_cov(end,end) = 20^2;
        % end
        b = thresholdTest(p.td_lambda,cov_minus, res_all, H_os, R);
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla');
        R_pva = [R_e2g, zeros(3,6);
            zeros(3,3), R_e2g, zeros(3,3);
            zeros(3,6), R_e2g];
        Rot_e2g = [R_pva, zeros(9,length(x_minus)-9);
            zeros(length(x_minus)-9, 9), eye(length(x_minus)-9)];
        H_os = H_os * Rot_e2g';
        cov_prior = Rot_e2g' * cov_minus * Rot_e2g;
        [x_plus,cov_plus,infor_ned,p.augcost] = ...
            mapUpdate(b, x_minus, cov_prior, res_all, H_os, R, Rot_e2g);
        p.infor_ned = infor_ned(1:length(x_plus)-length(lambda),1:length(x_plus)-length(lambda));
        p.num_meas_used = sum(b);
    case p.raps_ned_est
        % Solve in NED frame
        lla_deg = ecef2lla(x_minus(1:3)', 'WGS84');
        R_eg=computeRotForEcefToNed(lla_deg);
        if p.state_mode == p.pva_mode
            R_pva = [R_eg, zeros(3,6);
                zeros(3,3), R_eg, zeros(3,3);
                zeros(3,6), R_eg];
            [flag_rapid,p.num_sats_window] = checkRapidNumSatChange(p.num_sats_window, sum(cpt.num_sv~=0));
            % Rapid num of sat detected, may entering an open sky area,
            % reset prior covariance
            if flag_rapid == true
                p.state_cov = zeros(size(p.state_cov));
                p.state_cov(1:end-1,1:end-1) = diag(150^2*ones(1,length(x_minus)-1));
                p.state_cov(end,end) = 20^2;
            end
        else
            R_pva = R_eg;
        end
        Rot_e2g = [R_pva, zeros(num_user_states,length(x_minus)-num_user_states);
            zeros(length(x_minus)-num_user_states, num_user_states), eye(length(x_minus)-num_user_states)];
        Ht = H_os * Rot_e2g';
        xt_minus = Rot_e2g*(x_minus - x_minus);
        Pt_minus = Rot_e2g*p.state_cov*Rot_e2g';
        if p.state_mode == p.pva_mode
            num_constrain = 6;
            cov_spec_ecef = diag([p.raps.poshor_cov_spec; ...
                p.raps.poshor_cov_spec; p.raps.posver_cov_spec;...
                p.raps.velhor_cov_spec; p.raps.velhor_cov_spec;...
                p.raps.velver_cov_spec]);
            p_clk = diag([p.raps.va_cov_spec*ones(3,1);...
                p.raps.clk_cov_spec*ones(length(x_clk),1);...
                p.raps.dclk_cov_spec]);
            p_u = [cov_spec_ecef, zeros(6, length(x_minus)-6);
                zeros(length(x_minus)-6, 6), p_clk];
        elseif p.state_mode == p.pos_mode
            num_constrain = 3;
            cov_spec_ecef = diag([p.raps.poshor_cov_spec; ...
                p.raps.poshor_cov_spec; p.raps.posver_cov_spec]);
            p_clk = diag([p.raps.clk_cov_spec*ones(length(x_clk),1);...
                p.raps.dclk_cov_spec]);
            p_u = [cov_spec_ecef, zeros(3, length(x_minus)-3);
                zeros(length(x_minus)-3, 3), p_clk];
        end
        J_l = p_u^(-1);
        [flag,x_ned,cov_ned,b,J_out,p.augcost,num_nodes,constraint] = ...
            mapRiskAverse(num_constrain,res_all,Ht,Pt_minus,R,...
            diag(diag(J_l)),xt_minus,p.td_lambda,length(y_rho));
        cov_plus = Rot_e2g' * cov_ned * Rot_e2g;
        p.num_meas_used = sum(b);
        if p.state_mode == p.pva_mode
            p.raps_num_sat = sum(b(1:length(y_rho)));
        else
            p.raps_num_sat = sum(b);
        end
        p.infor_ned = J_out;
        p.raps_num_nodes = num_nodes;
        p.constraint = constraint;
        p.raps_flag = flag;
        x_plus = Rot_e2g'*x_ned + x_minus;
    otherwise
        error('Incorrect state estimation mode configuration');
end

p.state0 = x_plus(1:length(x_plus)-length(lambda));
p.state_cov = cov_plus(1:length(x_plus)-length(lambda),1:length(x_plus)-length(lambda));

% HH = diag(b)*H_os;
% HH([1,5,7],:) = [];
% HH = HH(:, 1:4);
% hSqrtInv = (HH'*HH)^(-1);
% PDOP = sqrt(trace(hSqrtInv(1:3,1:3)));
% if PDOP > 5
%     estState.pos = [];
%     return;
% end
res = res_all(1:num_rho);
estState.pos = x_plus(1:3);
if p.state_mode == p.pva_mode
    estState.vel = x_plus(4:6);
end
estState.clock_bias = x_plus(num_user_states+1);
estState.clock_drift = x_plus(end);

clk_est = x_plus(num_user_states+1:end-1);
j = 1;
for i = 1:length(cpt.num_sv)
    if cpt.num_sv(i) == 0
        continue;
    end
    estState.clock_sys(i) = clk_est(j);
    j=j+1;
end

end


