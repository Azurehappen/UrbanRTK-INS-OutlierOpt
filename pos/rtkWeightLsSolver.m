function [estState,res,state_cov] = rtkWeightLsSolver(p,cpt,code_only_flag)
% This is solver for computing user position, receiver clock
% bias, satellite system offset.
% Measurement selection applied
% Input: 
%       s_pos_ecef: 3-by-N Satellite position in ECEF frame.
%       x0 : 3-by-1 initial interative coordinates in ECEF frame.
%       y: m-by-1 Corrected pseudorange.
%
% Output:
%       
%       
%       
%       

%-------------------%
% Initialize
estState.clock_sys = dictionary;
estState.clock_sys(p.gps.sys_num) = NaN;
estState.clock_sys(p.glo.sys_num) = NaN;
estState.clock_sys(p.gal.sys_num) = NaN;
estState.clock_sys(p.bds.sys_num) = NaN;

x0 = p.state0(1:3);
[H_clk,x_clk] = formClkStatesAndH(cpt.num_sv);
xk = [x0;x_clk];

y_rho = cpt.corr_range;
num_rho = length(y_rho); % The number of measurement
y_phi = [];
num_phi = 0;
if p.post_mode == p.mode_rtkfloat && code_only_flag == false
    y_phi = cpt.phase_m;
    num_phi = length(y_phi);
    % if length(y_rho) ~= length(y_phi)
    %     warning('Num of pseudorange != phase');
    % end
end

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

if p.post_mode == p.mode_rtkfloat && code_only_flag == false
    xk = [xk;zeros(length(lambda),1)];
end

trans_mat = [eye(num_rho),zeros(num_rho,num_phi);
            zeros(size(T_dd,1),num_rho),T_dd];

H_rho = zeros(num_rho, length(xk));
H_rho(:, 4:(4+length(x_clk)-1)) = H_clk;
H_phi = zeros(num_phi, length(xk));

s_pos_ecef = cpt.s_pos_ecef;
if p.post_mode == 1 && p.IGS_enable == 1
    s_pos_ecef = cpt.s_pos_prc;
end

Range = zeros(num_rho,1);
geo_range_rho = zeros(num_rho,1);
geo_range_phi = zeros(num_phi,1);
meas = [y_rho;T_dd * y_phi];
if code_only_flag == true
    W = eye(num_rho);
else
    [R_rho, R_phi] = constructMeasNoise(p, cpt, 1);
    W_rho = diag(1./diag(R_rho));
    Rt = T_dd * R_phi * T_dd';
    W_phi = Rt^(-1);
    W = [W_rho,zeros(size(W_rho,1),size(W_phi,2));
        zeros(size(W_phi,1),size(W_rho,2)),W_phi];
end
for iter=1:p.Nls
    for j=1:num_rho
        Range(j)=norm(s_pos_ecef(:,j)-xk(1:3));
        los = (xk(1:3)-s_pos_ecef(:,j))'/Range(j)+...
            [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0];
        H_rho(j,1:3)=los;
        geo_range_rho(j) = Range(j)+sagnac(p,s_pos_ecef(:,j),xk);
        if p.post_mode == p.mode_rtkfloat && code_only_flag == false
            H_phi(j,1:3)=los;
            geo_range_phi(j) = geo_range_rho(j);
        end
    end
    if p.post_mode == p.mode_rtkfloat && code_only_flag == false
        geo_range_phi = T_dd * geo_range_phi;
        H_phi_dd = T_dd * H_phi;
        H_phi_dd(:, length(xk)-length(lambda)+1:end) = diag(lambda);
        H = [H_rho;H_phi_dd];
        geo_range = [geo_range_rho;geo_range_phi];
    else
        H = H_rho;
        geo_range = geo_range_rho;
    end
    res = meas - geo_range - H(:, 4:end) * xk(4:end);
    delta_x = (H'*W*H)^(-1)*H'*W*(res);
    xk=xk+delta_x;
    if (norm(delta_x(1:4)) < p.LSthrsh)
        break;
    end
    if (iter>p.Nls)&& (norm(delta_x) > p.LSthrsh)
        warning('Postion path length iteration failed in user_pos calculation');
    end
end
%------------------%
% [estState.pos,estState.clock_bias,isb_est,res] = LSsolver(p,xk,H_isb,cpt);
estState.pos = xk(1:3);
estState.clock_bias = xk(4);
clk_est = xk(4:4+length(x_clk)-1);

state_cov = (H'*W*H)^(-1);
res = res(1:num_rho);
j = 1;
for i = 1:length(cpt.num_sv)
    if cpt.num_sv(i) == 0
        continue;
    end
    estState.clock_sys(i) = clk_est(j);
    j=j+1;
end

end
