function [res_dd, H_dd, R_dd] = doubleDiff(cpt, res_all, H_os, R, min_sv)

dd_T = [];
row_ind = 0;
to_delete = [];
num_sv = cpt.num_sv;
num_sv_re = [];
for i = 1:length(num_sv)
    if num_sv(i) == 0
        continue;
    end
    %Find the pivot satellite
    ind_range = row_ind+1:row_ind+num_sv(i);
    elev = cpt.elev(ind_range);
    pivot = find(elev== max(elev));
    row_ind = row_ind + num_sv(i);
    if elev(pivot) < deg2rad(55)
        to_delete = [to_delete,ind_range];
        continue;
    end
    len = num_sv(i);
    num_sv_re = [num_sv_re,len-1];
    dd_sub = eye(len);
    dd_sub(:,pivot) = -1;
    dd_sub(pivot,:) = [];
    dd_T = blkdiag(dd_T,dd_sub);
end

if isempty(dd_T) || sum(num_sv_re)<min_sv
    % No positioning will be performed
    res_dd = [];
    H_dd = [];
    R_dd = [];
end

dd_T = blkdiag(dd_T,dd_T);
meas_l = length(res_all)/2;
to_delete = [to_delete,meas_l+to_delete];
res_all(to_delete) = [];
H_os(to_delete,:) = [];
R(to_delete,:) = [];
R(:,to_delete) = [];

res_dd = dd_T*res_all;
H_dd = dd_T*H_os;
R_dd = dd_T*R*dd_T';

end