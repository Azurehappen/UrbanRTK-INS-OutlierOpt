function cpt = cleanMeasDataStruct(cpt, del_ind)

% Get the field names
fields = fieldnames(cpt);

% Loop over each field
for i = 1:length(fields)
    % Get the field name
    field = fields{i};

    if strcmp(field,'num_sv') || strcmp(field,'gps_range') ||...
       strcmp(field,'gps_sat_pos') || strcmp(field,'is_open_sky')
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