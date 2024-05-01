function tidx = ephtidx(eph_info,t_sv,prn,message_duration,sys_type,ppp_flag)
% find the time index in eph data

t_oc = eph_info.t_oc{prn};
SV_health = eph_info.SV_health(prn,:);
dtr = t_sv-t_oc;
if ppp_flag == true
    tidx = find(dtr>=-message_duration&dtr<=message_duration);
else
    tidx = find(dtr>=-message_duration&dtr<=message_duration);
end
% Satellite health check
i = SV_health(tidx)==0;
tidx = tidx(i);
if strcmp(sys_type,'gal')
    % Skip F/NAV msg
    data_source = eph_info.Data_source(prn,:);
    % For Septentrio receivers, use F/NAV msg (258), otherwise 517
    i = data_source(tidx)==258;
    tidx = tidx(i);
end
end