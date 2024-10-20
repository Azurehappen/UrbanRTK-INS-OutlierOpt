function [x_plus,error_x_plus,cov_plus,J_ned,aug_cost] = mapUpdate(p,by, ...
    x_prior, err_x_prior, cov_prior, meas_res, H, R, Rot_e2g)
% Maximum A Posteriori state estimate and the augmented cost function
% when doing measurement selection
% INPUT: by  - measurement selection vector
%        E_P - sqrt(P^-1)
%        E_R - sqrt(R^-1)

E_P = chol(cov_prior^-1);    % cholesky decomp. of state prior info. matrix
E_R = chol(R^-1);    % cholesky decomp. of measurement info. matrix

if isempty(by)
    error('Input arg "by" is empty.');
else
    Phiby = diag(by);
    A = [E_R * Phiby * H; E_P];             
    c = [E_R * Phiby * meas_res; zeros(length(err_x_prior),1)];   
    dx = (A'*A)^-1*A'*c; % LeastSquares % posterior state estimate
    aug_cost = norm(A*dx-c)^2;
end

if p.state_mode == p.ins_mode
    x_plus = x_prior;
    % correct pos and vel states
    x_plus(1:6) = x_prior(1:6) + Rot_e2g(1:6,1:6)' * dx(1:6);
    % Correct the quaternion (Farrell Appendix.D Exercise D.4)
    rho = dx(7:9);
    rho_bar = 0.5*rho;
    % q_u= sqrt(1-(norm(rho_bar))^2);
    % q_e=[q_u;rho_bar];
    if norm(rho_bar) > 1
        q_e = 1/sqrt(1+norm(rho_bar)^2)*[1; rho_bar];
        % warning('norm(rho_bar) > 1');
    else
        q_e = [sqrt(1-norm(rho_bar)^2); rho_bar];
    end
    x_plus(7:10) = quatInv(quatMult(q_e,quatInv(x_prior(7:10,1))));
    % Correct other states
    x_plus(11:end) = x_prior(11:end) + dx(10:end);
else
    x_plus(1:9) = x_prior(1:9) + Rot_e2g(1:9,1:9)'*dx(1:9);
    x_plus(10:end) = x_prior(10:end) + dx(10:end);
end
error_x_plus = dx;
PhiH = Phiby * H;
J_ned    = PhiH' * R^-1 * PhiH + cov_prior^-1; % eqn. (13) in [1]
cov_ned = J_ned^-1;
cov_plus = Rot_e2g' * cov_ned * Rot_e2g;

end