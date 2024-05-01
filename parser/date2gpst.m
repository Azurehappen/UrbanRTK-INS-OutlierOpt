function [gps_week, gps_dow, gps_sow] = date2gpst(date)
% Convert from calendar date to gnss system time
% The number of day of the starting time from each system to GPS is an
% integer multiple of 7, it use the gps time week num as reference.
% Parameters: 
%           date: [year, month, day, hour, min, sec]
%
% Return Values:
%   week:   Week number
%   sow:    seconds of week
%   dow:    day of week
%
% Author: Azurehappen

% Getting starting time date number
start_datenum = 723186; % GPS start Date number: datenum([1980,1,6,0,0,0])
%number of days since the beginning of GPS time
deltat   = (datenum([date(:,1), date(:,2), date(:,3)]) - start_datenum);
gps_week = floor(deltat/7);
gps_dow  = deltat - gps_week*7;
gps_sow  = (deltat - gps_week*7)*86400; 
gps_sow = gps_sow + date(:,4)*3600 + date(:,5)*60 + date(:,6);

        