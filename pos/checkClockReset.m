function [x, cov, flag] = checkClockReset(p, x0, err_cov, num_user_states, res, cpt)

x = x0;
cov = err_cov;
flag = false;
% err_cov is the error state covariance
% x0 is the state. len(x0)-len(err_state) = 1 because 4 Quat compare to 3
% attitude error states
if sum(abs(res) > 5000) > 0.8*length(res)
    [estState,~] = weightLsSolver(p,cpt,true);
    if max(cpt.num_sv(cpt.num_sv>0)) <= 4%p.min_sv
        warning('No. of sat too less, cannot estimate the clock bias');
    end
    x(num_user_states+1:end-1) = estState.clock_bias;
    x(end) = estState.clock_drift;
    cov(:, num_user_states-1+1:end) = zeros(length(x0)-1,length(x0)-num_user_states);
    cov(num_user_states-1+1:end,:) = zeros(length(x0)-num_user_states,length(x0)-1);
    clk_len = length(x0) - num_user_states - 1;
    cov(num_user_states-1+1:end-1, num_user_states-1+1:end-1) = 200^2*diag(ones(1,clk_len));
    cov(end, end) = 5^2;
    flag = true;
end

end