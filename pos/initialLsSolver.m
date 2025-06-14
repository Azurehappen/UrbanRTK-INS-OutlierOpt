function [estState,res] = initialLsSolver(p,cpt)
% This is solver for computing user position, receiver clock
% bias, satellite system offset.
% Measurement selection applied
% Input: 
%       s_pos_ecef: 3-by-N Satellite position in ECEF frame.
%       x0 : 3-by-1 initial iterative coordinates in ECEF frame.
%       y: m-by-1 Corrected pseudorange.
%
% Output:

%-------------------%
% Initialize
estState.clock_sys = dictionary;
estState.clock_sys(p.gps.sys_num) = NaN;
estState.clock_sys(p.glo.sys_num) = NaN;
estState.clock_sys(p.gal.sys_num) = NaN;
estState.clock_sys(p.bds.sys_num) = NaN;

[H_clk,x_clk] = formClkStatesAndH(cpt.num_sv);
xk = [zeros(3,1);x_clk]; % POS, CLK

y_rho = cpt.corr_range;
num_rho = length(y_rho); % The number of measurement

H_rho = zeros(num_rho, length(xk));
H_rho(:, 4:(4+length(x_clk)-1)) = H_clk;

s_pos_ecef = cpt.s_pos_ecef;
if p.post_mode == p.mode_ppp && p.IGS_enable == 1
    s_pos_ecef = cpt.s_pos_prc;
end

Range = zeros(num_rho,1);
geo_range_rho = zeros(num_rho,1);
meas = y_rho;
for iter=1:p.Nls
    for j=1:num_rho
        Range(j)=norm(s_pos_ecef(:,j)-xk(1:3));
        los = (xk(1:3)-s_pos_ecef(:,j))'/Range(j)+...
            [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0];
        H_rho(j,1:3)=los;
        geo_range_rho(j) = Range(j)+sagnac(p,s_pos_ecef(:,j),xk);
    end
    H = H_rho;
    geo_range = geo_range_rho;
    res = meas - geo_range - H(:, 4:end) * xk(4:end);
    delta_x = (H'*H)^(-1)*H'*(res);
    xk=xk+delta_x;
    if (norm(delta_x) < p.LSthrsh)
        break;
    end
    if (iter>p.Nls)&& (norm(delta_x) > p.LSthrsh)
        warning('Postion path length iteration failed in user_pos calculation');
    end
end
%------------------%
estState.pos = xk(1:3);
estState.clock_bias = xk(4);
clk_est = xk(4:4+length(x_clk)-1);

estState.state_cov = (H'*H)^(-1);
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
