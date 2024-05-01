function ECEF_del = RAC2ECEF(RAC_del,orbit_pos,orbit_vel)
%%% Transform the coordinate from Radial,Along-track,Cross-track frame to
%%% ECEF frame.
    %%% RAC_del is the RAC coordiante 3*1
    %%% orbit_pos is the orbit_position in ECEF of the SV 3*1
    %%% orbit_vel is the velocity in ECEF of the SV on the orbit 3*1
   
    %%%%% eqn.2(Hadas-Bosy2015)
    unit_a = orbit_vel/norm(orbit_vel);%unit vector in radial
    unit_c = cross(orbit_pos,orbit_vel)/norm(cross(orbit_pos,orbit_vel));%unit vector in along-track
    unit_r = cross(unit_a,unit_c);%unit vector in cross-track
    
    %%%%% eqn.3(Hadas-Bosy2015)
    ECEF_del = [unit_r,unit_a,unit_c]*RAC_del;
    
    


end

%%Reference:
% IGS RTS precise orbits and clocks evrification and quality degradation
% over time (GPS solutions Jan 2014)