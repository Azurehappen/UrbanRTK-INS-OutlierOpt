function [x_plus,error_x_plus,cov_plus,J_info,aug_cost] = ekfUpdate(x_prior, err_x_prior,cov_prior, meas_res, H, R)

E_P = chol(cov_prior^-1);    % cholesky decomp. of state prior info. matrix
E_R = chol(R^-1);    % cholesky decomp. of measurement info. matrix

A = [E_R * H; E_P];             
c = [E_R * meas_res; zeros(length(err_x_prior),1)];   
dx_map = (A'*A)^-1*A'*c; % LeastSquares % posterior state estimate
aug_cost = norm(A*dx_map-c)^2;

% Innovation (or residual) covariance
S = H * cov_prior * H' + R;
% Near-optimal Kalman Gain
Kk = cov_prior * H' * S^(-1);
dx = Kk * meas_res;
I = eye(size(cov_prior));
cov_plus = (I - Kk*H)*cov_prior*(I - Kk*H)' + Kk*R*Kk';

x_plus = x_prior;
% correct pos and vel states
x_plus(1:6) = x_prior(1:6) + dx(1:6);
% Correct the quaternion (Farrell Appendix.D Exercise D.4)
rho = dx(7:9);
rho_bar = 0.5*rho;
% q_u= sqrt(1-(norm(rho_bar))^2);
% q_e=[q_u;rho_bar];
if norm(rho_bar) > 1
    q_e = 1/sqrt(1+norm(rho_bar)^2)*[1; rho_bar];
    warning("Large ")
else
    q_e = [sqrt(1-norm(rho_bar)^2); rho_bar];
end
x_plus(7:10) = quatInv(quatMult(q_e,quatInv(x_prior(7:10,1))));
% Correct other states
x_plus(11:end) = x_prior(11:end) + dx(10:end);
error_x_plus = dx;
J_info    = cov_plus^(-1); % eqn. (13) in [1]