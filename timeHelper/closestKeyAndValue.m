function [closestKey, val] = closestKeyAndValue(dict, x, low_bound, upper_bound)
    keys = dict.keys; % Convert keys to array of doubles
    diffs = x - keys; % Compute differences
    valid_indices = find(diffs >= low_bound & diffs <= upper_bound);
    
    if isempty(valid_indices) % If no diffs are within bound
        closestKey = [];
        val = [];
    else
        [~, idx] = min(diffs(valid_indices)); % Find minimum diff index within valid indices
        closestKey = keys(valid_indices(idx));
        val = dict(closestKey);
    end
end