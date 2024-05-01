function cpt = vrs_corr_compute(p,cpt,eph,tdoy,rt)

[p.lat, p.lon, p.h_r, ~, ~, ~] = ecef2llh(p,p.P_base);
load geoid;
%%%%% geodetic ellipsoidal separation
N = ltln2val(geoid, geoidrefvec, p.lat*180/pi, p.lon*180/pi);
%%%%% orthometric height of the receiver (Groves 2.121)
H_r=p.h_r-N;
vrs_psedR = NaN(length(cpt.corr_range),1);
diff_corr = NaN(length(cpt.corr_range),1);
sys_i = find(cpt.svprn_mark~=0);
obs = [];
for i = 1:length(cpt.corr_range)
    sys_num = cpt.svprn_mark(sys_i(i));
    switch sys_num
        case 1
            eph_info = eph.GPS;
            sys = 'GPS';
        case 2
            eph_info = eph.GLO;
            sys = 'GLO';
            rt.sow = time_shift(rt.sow - p.glo.lps_gps);
        case 3
            eph_info = eph.GAL;
            sys = 'GAL';
        case 4
            eph_info = eph.BDS;
            sys  = 'BDS';
            rt.sow = time_shift(rt.sow - p.bds.lps_gps);
    end
    prn = cpt.prn_record(sys_i(i));
    tidx = ephtidx(eph_info.t_oc{prn},rt.sow,eph_info.SV_health(prn,:),p.gps.message_duration);
    p.post_mode = 1;p.IGS_enable = 1;
    if ~isempty(tidx)
        [tp,dt_sv,sat] = inverse_eph(p,eph_info,obs,prn,tidx,rt.sow,sys,p.P_base);
        if ~isnan(tp)
%             theta = p.omge * (tp+dt_sv);
%             R = [ cos(theta)  sin(theta)  0;   
%                 -sin(theta)  cos(theta)  0; 
%                 0               0       1];
%             sat_pos_Rcorr = R*sat.pos_prc;
%             sat_ecef_Rcorr = R*sat.pos_ecef;
            [trop_delay, ~, ~, ~, ~]=UNB3M(p.lat,H_r,tdoy,cpt.elev(i));
            [iono_delay] = ustec_iono_delay_computation(p,p.USTEC,cpt.elev(i),cpt.az(i),rt,p.L1freq);
            vrs_psedR(i) = norm(p.P_base-sat.pos_prc)+sagnac(p,sat.pos_prc,p.P_base)-dt_sv*p.c+trop_delay+iono_delay;
        else
            vrs_psedR(i) = NaN;
            cpt.num_sv(sys_num) = cpt.num_sv(sys_num)-1;
        end        
    end
    if ~isnan(vrs_psedR(i))
        p.post_mode = 0;p.IGS_enable = 0;
        tidx = ephtidx(eph_info.t_oc{prn},rt.sow,eph_info.SV_health(prn,:),p.gps.message_duration);
        t_sv = rt.sow - vrs_psedR(i)/p.c;
        [sat, dt_sv] = eph2pos(p,eph_info,obs,prn,tidx,t_sv,sys);
        range = norm(p.P_base-sat.pos_ecef)+sagnac(p,sat.pos_ecef,p.P_base);
        diff_corr(i) = vrs_psedR(i) + p.c*dt_sv - range;
    end
end
% Delete the data that has no diff correction
del_ind = find(isnan(diff_corr));

if ~isempty(del_ind)
    ind_prn = find(cpt.prn_record~=0);
    cpt.prn_record(ind_prn(del_ind))=0;
    cpt.svprn_mark(ind_prn(del_ind))=0;
    diff_corr(del_ind) = [];
    cpt.corr_range(del_ind) = [];
    cpt.s_pos_ecef(:,del_ind) = [];
    cpt.s_v_ecef(:,del_ind) = [];
    cpt.tp(del_ind) = [];
    cpt.elev(del_ind) = [];
    cpt.az(del_ind) = [];
    cpt.sat_pos_Rcorr(:,del_ind) = [];
    cpt.sat_posprc_Rcorr(:,del_ind) = [];
    cpt.sat_v_Rcorr(:,del_ind) = [];
end
cpt.diff_corr = diff_corr;

end