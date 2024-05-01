function clock_corr_m = obtainSatClockCorr(clock_dict, t, low_bound, ...
    upper_bound, prn, iod, sys_type)
% Obtain precise satellite clock bias
keys = clock_dict.keys; % Convert keys to array of doubles
diffs = t - keys; % Compute differences
valid_indices = find(diffs >= low_bound & diffs <= upper_bound);
clock_corr_m = NaN;
if isempty(valid_indices)
    return;
end

[~, idx] = min(diffs(valid_indices)); % Find minimum diff index within valid indices

while idx > 0
    key = keys(valid_indices(idx));
    val = clock_dict(key);
    sat_clock = val.(sys_type)(prn);
    dt = t - key;
    if sat_clock.iod == iod
        clock_corr_m = sat_clock.corr + sat_clock.vel*dt + sat_clock.acc*dt^2;
        return;
    end
    idx = idx - 1;
end

end