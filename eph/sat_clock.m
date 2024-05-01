function dt_sv = sat_clock(sysp, prn, tidx, eph, Ek, t_sv)
% Compute satellite clock correction (m) for pseudorange processing
% Reference: 
% 'https://www.gps.gov/technical/icwg/IS-GPS-200H.pdf'
% 'https://www.gsc-europa.eu/sites/default/files/sites/all/files/Galileo-OS-SIS-ICD.pdf'
% 'http://en.beidou.gov.cn/SYSTEMS/ICD/201902/P020190227702348791891.pdf'
% where, F = -2*sqrt(mu)/c^2
% SV PRN code phase time offset (s), account for week rollovers
dt_oc = limit_tgps(t_sv - eph.t_oc{prn}(tidx));
    
% relativistic correction (s)
% dt_R = sysp.F*eph.e(prn,tidx)*eph.sqrtA(prn,tidx)*sin(Ek);

% SV PRN code phase offset (s)
dt_sv = eph.a_f0(prn,tidx) + eph.a_f1(prn,tidx)*dt_oc + eph.a_f2(prn,tidx)*dt_oc^2;% + dt_R;

end