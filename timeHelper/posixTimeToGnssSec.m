function gnss_sow = posixTimeToGnssSec(t)

a = datetime(t, 'Convertfrom', 'posixtime');
datet_array = [year(a),month(a),day(a),...
    hour(a),minute(a),second(a)];
[~, ~, gnss_sow] = date2gnsst(datet_array);

end
