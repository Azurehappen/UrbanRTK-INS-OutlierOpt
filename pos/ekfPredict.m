function [prior_state, prior_cov] = ekfPredict(p, state, cov, dt)

num_user_states = p.modeToNumUserErrStates(p.state_mode);
num_clk = length(state) - num_user_states - 1;
Qc = p.ekf_para.q_clkDrift * [dt^3/3*ones(num_clk,num_clk), dt^2/2*ones(num_clk,1);
                   dt^2/2*ones(1,num_clk), dt];
Qp_diag = p.ekf_para.q_pos * ones(1,3);
Qv_diag = p.ekf_para.q_vel * ones(1,3);

lla_deg = ecef2lla(state(1:3)', 'WGS84');
Rot_eg=computeRotForEcefToNed(lla_deg);
Qa_ned = diag([p.ekf_para.q_accHor,p.ekf_para.q_accHor,p.ekf_para.q_accVer]);
Qa_ecef = Rot_eg*Qa_ned*Rot_eg';
phi = diag(ones(1, length(state)));
if p.state_mode == p.pva_mode
    % x = [x,y,z,vx,vy,vz,ax,ay,az,clks,clk_drift]
    phi(1:9,1:9) = [eye(3,3),dt*eye(3,3), dt^2/2 * eye(3,3);
                zeros(3,3),eye(3,3),dt*eye(3,3);
                zeros(3,3),zeros(3,3),eye(3,3)];
    % see Aided Nav. Book eqn. 4.110
    Qd = [Qa_ecef*dt^5/20, Qa_ecef*dt^4/8, Qa_ecef*dt^3/6;
           Qa_ecef*dt^4/8, Qa_ecef*dt^3/3, Qa_ecef*dt^2/2;
           Qa_ecef*dt^3/6, Qa_ecef*dt^2/2, Qa_ecef*dt];
elseif p.state_mode == p.pos_mode
     % x = [x,y,z,clks,clk_drift]
    Qd = diag(Qp_diag);
end
if p.double_diff == false
    Qd = [Qd, zeros(num_user_states, num_clk+1);
        zeros(num_clk+1, num_user_states), Qc];
    phi(num_user_states+1:length(state)-1, end) = dt; % clk(k) = clk(k-1) + clk_drift * dt
end


prior_state = phi * state;
prior_cov = phi * cov * phi' + Qd;