function [flag_pos,flag_vel] = rapsValidation(td_lambda,cov_prior, meas_res, H, R, b_raps, num_y)
% Check the reasidual threshold to validate the RAPS results
% If in the extreme case that most of measurements have outliers,
% RAPS would have to select those measurements to satisfy the constraint.
% This will only be used to PVA model
b = thresholdTest(td_lambda,cov_prior, meas_res, H, R);
b_range = b(1:num_y); % Check range measurements only
b_raps_range = b_raps(1:num_y);
br_selcted = b_range(b_raps_range==1);
if sum(b_range) < 0.8*length(b_range) || ~all(br_selcted == 1)
    flag_pos = false;
else
    flag_pos = true;
end
b_vel = b(num_y+1:end); % Check range measurements only
b_raps_vel = b_raps(num_y+1:end);
bv_selcted = b_vel(b_raps_vel==1);
if sum(b_vel) < 0.8*length(b_vel) || ~all(bv_selcted == 1)
    flag_vel = false;
else
    flag_vel = true;
end

end