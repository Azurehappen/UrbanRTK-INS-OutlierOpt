function cpt = elevaz_check(p,cpt,u_pos)
% computing the azimuth and elevation angle from user to satellite (s)
% Remove the satellite values that not satisfy elevation require
%
% Input
%       p: paramters
%       cpt: computed data of svprn_mar k, corr_range, s_pos_ecef,num_sv
%       u_pos estimated user position
ind_prn = find(cpt.svprn_mark~=0);
num = length(cpt.num_sv);
count = 0;
gps_range = [];
gps_sat_pos = [];
for i = 1:num
    len = cpt.num_sv(i);
    for j = 1:len
        % Compute the elev and az
        [cpt.elev(count+j), cpt.az(count+j)] = sat_elev_azimuth(p,u_pos,cpt.sat_pos_Rcorr(:,count+j));
        if cpt.elev(count+j) >= p.open_sky_elev_rad && cpt.svprn_mark(count+j) == p.gps.sys_num
             gps_range = [gps_range;cpt.corr_range(count+j)];
             gps_sat_pos = [gps_sat_pos,cpt.sat_pos_Rcorr(:,count+j)];
        end
        if cpt.elev(count+j) < p.elev_mark_rad
            % Eliminate the variables where not yield the elev mask,
            cpt.corr_range(count+j) = 0;
            cpt.doppler(count+j) = 0;
            cpt.svprn_mark(ind_prn(count+j)) = 0;
            cpt.prn_record(ind_prn(count+j)) = 0;
            cpt.num_sv(i) = cpt.num_sv(i)-1;
        end
    end
    count = count + len;
end

% Get the field names
fields = fieldnames(cpt);

% Find indices to delete
del_ind = find(cpt.corr_range==0);

% Loop over each field
for i = 1:length(fields)
    % Get the field name
    field = fields{i};

    if strcmp(field,'num_sv')
        continue;
    end

    % Check the size of the field
    if size(cpt.(field), 1) == 3
        % If the field is pos, vel...
        cpt.(field)(:, del_ind) = [];
    else
        cpt.(field)(del_ind, :) = [];
    end
end
cpt.gps_range = gps_range;
cpt.gps_sat_pos = gps_sat_pos;

end