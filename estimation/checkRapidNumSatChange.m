function [flag,num_sat_array] = checkRapidNumSatChange(num_sat_array,num_sys)

max_num = max(num_sat_array);
min_num = min(num_sat_array);
if max_num - min_num >= num_sys*2 && num_sat_array(end) == max_num
    flag = true;
    num_sat_array(:) = NaN;
    return;
end
flag = false;

end