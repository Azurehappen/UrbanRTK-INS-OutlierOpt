function b = thresholdTest(td_lambda,cov_prior, meas_res, H, R)

% Compute the covariance of residual component
sigma_r = zeros(length(meas_res),1);
for i=1:length(meas_res)
    sigma_r(i) =  H(i,:)*cov_prior*H(i,:)' + R(i,i);
end
sigma_r = td_lambda * sigma_r;
b = meas_res < sigma_r;