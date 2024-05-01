function [sat, dt_sv, ddt_sv] = geph2pos(p,eph,prn,tidx,t_gst)
% compute ephemeris to satellite position and clock bias
% for GLONASS system
%%%%%% Inputs
% p - parameters
% prn - svid
% eph  - ephemeris data structure
% obs  - oberservable data structure
% prn  - satellite prn
% tidx - time index in eph(row)
% t_gst - transmit time
%%%%%% Outputs
% sat_pos_ecef - satellite position in ECEF
% sat_v_ecef - satellite velocity in ECEF
% dt_sv - satellite clock bias
% Reference: GLONASS ICD 5.1.5
sat.pos_ecef = NaN(3,1); 
sat.v_ecef = NaN(3,1);
dt_sv = 0;
dt_sv_p = 0;
if p.post_mode == p.mode_ppp
    % Save the precise pos and velocity in PPP mode 
    sat.pos_prc = NaN(3,1); 
%     sat.v_prc = NaN(3,1);
end
%Psedotime that the satellite broadcast the signal
%--------------------------%
% Find if there is PPP IGS correction
% [tidx,obt_idx,clk_idx,pos_tage,IGSdata,icb] = tidxconf(p,t_gst,prn,tidx,eph,sys_type);
tidx=tidx(end);
%--------------------------%
% position
X_pos = eph.X(prn,tidx);
Y_pos = eph.Y(prn,tidx);
Z_pos = eph.Z(prn,tidx);
% velocity
Xdot  = eph.Xdot(prn,tidx);
Ydot  = eph.Ydot(prn,tidx);
Zdot  = eph.Zdot(prn,tidx);
% acceleration 
Xacc  = eph.Xacc(prn,tidx);
Yacc  = eph.Yacc(prn,tidx);
Zacc  = eph.Zacc(prn,tidx);

Xo = [X_pos,Y_pos,Z_pos,Xdot,Ydot,Zdot];
acc = [Xacc,Yacc,Zacc];
% compute clock correction
dt_sv = eph.nTauN(prn,tidx)+eph.pGammaN(prn,tidx)*(t_gst-eph.t_oc{prn}(tidx));
time = t_gst-eph.t_oc{prn}(tidx)-dt_sv;
Xf = ode4(@(t,y) diffEq(t,y, acc), 0,60,time, Xo);
sat.pos_ecef = Xf(end,1:3)';
sat.v_ecef = Xf(end,4:6)';

dt_sv2 = eph.nTauN(prn,tidx)+eph.pGammaN(prn,tidx)*(t_gst+0.001-eph.t_oc{prn}(tidx));
ddt_sv = (dt_sv2 - dt_sv) / 0.001;

if p.post_mode == p.mode_ppp && p.IGS_enable == 1
    sat.pos_prc = sat_position_precise(p,IGSdata,sat.pos_ecef,sat.v_ecef,prn,obt_idx,t_gst);
    dt_sv_p = sat_clock_precise(p,IGSdata,prn,clk_idx,t_gst);
    dt_sv = dt_sv + dt_sv_p + icb(prn)*1e-9;
end

%---------------------------------------------------------%
function y_out = ode4(F,t0,h,tfinal,y0)
% ODE4  Classical Runge-Kutta ODE solver.
% t0  - initial time
% h  - time step
% tfinal  - simulation duration
%   to solve
%      dy/dt = F(t,y)
%   with y(t0) = y0.
    tt = t0;
    y_out = y0;
    if tfinal<0
        h = -60;
    end
    N = floor((tfinal-t0)/h);
    r = tfinal-t0-N*h;
    if N>0
        for i = 1:N
            s1 = F(tt+h,y_out);
            s2 = F(tt+h/2, y_out+h*s1/2);
            s3 = F(tt+h/2, y_out+h*s2/2);
            s4 = F(tt+h, y_out+h*s3);
            y_out = y_out + h*(s1 + 2*s2 + 2*s3 + s4)/6;
            tt = tt + h;
        end
    end
    if abs(r)>0
        s1 = F(tt+r,y_out);
        s2 = F(tt+r/2, y_out+r*s1/2);
        s3 = F(tt+r/2, y_out+r*s2/2);
        s4 = F(tt+r, y_out+r*s3);
        y_out = y_out + r*(s1 + 2*s2 + 2*s3 + s4)/6;
    end
end
%---------------------------------------------------------%
function dydt = diffEq(t,y, acc)
% Satellite positions integration equations. 
%
%   Syntax:
%       dy = diffEq (t, y, acc_sun, acc_moon)
%   Parameters:
%       y:     a 6-element column vector of positions and velocities          
%               [y]=[xo,yo,zo,Vxo,Vyo,Vzo]
%
%       acc_sun:    accelerations induced by gravitational perturbations of the Sun. 
%                   It is a 3-element vector. 
%       acc_moon:   accelerations induced by gravitational perturbations of the Moon. 
%                   It is a 3-element vector.        
%       p: parameters
%   Return values: 
%       dy:   The respective equations as specified in the ICD

mu = p.glo.mu;
a_e = p.glo.a_e;
C_20 = p.glo.C_20;
omega = p.glo.OmegaDot_e;

dydt = zeros(1,6);
dydt(1) = y(4);
dydt(2) = y(5);
dydt(3) = y(6);
r=norm(y(1:3)); 

dydt(4) = -mu*y(1)/(r^3) - 1.5*C_20*mu*(a_e^2/r^5)*y(1)*(1-5*y(3)^2/r^2) + omega^2*y(1) + acc(1) ...
        +2*omega*y(5);
dydt(5) = -mu*y(2)/(r^3) - 1.5*C_20*mu*(a_e^2/r^5)*y(2)*(1-5*y(3)^2/r^2) + omega^2*y(2) + acc(2) ...
        -2*omega*y(4);
dydt(6) = -mu*y(3)/(r^3) - 1.5*C_20*mu*(a_e^2/r^5)*y(3)*(3-5*y(3)^2/r^2) + acc(3);
end
%---------------------------------------------------------%

end