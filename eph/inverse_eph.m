function [tp,dt_sv,sat] = inverse_eph(p,eph,obs,prn,tidx,obs_tr,sys,re_pos)
% p: parameters
% tidx: column index for eph data
% obs_tr: reveiver time
% sys: system type
% re_pos: user position

% Initialize tm (time of propagation)
tp = 2.5e+07/p.c;
for i = 1:40
    [sat, dt_sv] = eph2pos(p,eph,obs,prn,tidx,obs_tr-tp,sys);
    if ~isnan(sat.pos_ecef(1))
    sat_pos = sat.pos_ecef; sat_v = sat.v_ecef;
    h = norm(re_pos - sat_pos)+sagnac(p,sat_pos,re_pos)-(tp+dt_sv)*p.c;
%     ddt_sv = -eph.a_f1(prn,tidx)-2*eph.a_f2(prn,tidx)*limit_tgps(obs_tr-tp- eph.t_oc(prn,tidx));
%     dh = -sum(sat_v)/norm(re_pos - sat_pos)-...
%         p.omge*(sat_v(1)*re_pos(2)-sat_v(2)*re_pos(1))/p.c-p.c-ddt_sv*p.c; 
    %'ddt_sv*p.c' is very small, don't have to include this, but it's ok to keep it.
    dh = -sum(sat_v.*(sat_pos-re_pos))/norm(re_pos - sat_pos)+...
        p.omge*(sat_v(1)*re_pos(2)-sat_v(2)*re_pos(1))/p.c-p.c;
    tm_old = tp;
    tp = tp - h/dh;
    if abs(tp-tm_old)<10e-11
        break;
    end
    else
        tp = NaN;
        break;
    end
end

end