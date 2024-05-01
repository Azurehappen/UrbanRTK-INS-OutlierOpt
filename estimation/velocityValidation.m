function flag = velocityValidation(prior_vel_ned,vel_est_ned,prior_ver_p)
% Check the reasidual threshold to validate the velocity
% The vertical velocity shouldn't have weird magnitude (> 5m/s)
if abs(vel_est_ned(3)) > 5 %|| abs(prior_vel_ned(3) - vel_est_ned(3)) > 2*sqrt(prior_ver_p)
    flag = false;
    return;
end
flag = true;

end