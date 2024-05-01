function [H_isb,x_isb] = formIsbStatesAndH(p, num_sys, mode)
% Create offset matrix for multi-GNSS
% Input:
%       num_offset: 4-1 vector, the number of satellite for each system
%       num_sys: 
% Outout:
%       offset: partial H matrix for offset paramters
total = sum(num_sys);
isb_status = [p.enableGLO, p.enableGAL, p.enableBDS];
if mode == p.mode_sps
    num_offset = sum(num_sys~=0)-1;
    ind = find(num_sys~=0);
    H_isb = zeros(total,num_offset);
    start = 0;
    if ~isempty(H_isb) % If only one system, don't need offset
        for i = 1:num_offset
            start = start + num_sys(ind(i));
            H_isb(start+1:start+num_sys(ind(i+1)),i)=1;
        end
    end
    x_isb = zeros(num_offset,1); % offset parts for state vector
else
    num_isb = p.enableGLO + p.enableGAL + p.enableBDS;
    x_isb = zeros(num_isb,1);
    H_isb = zeros(total,sum(num_isb));
    if ~isempty(H_isb)
        for i = 1:length(isb_status)
            if isb_status(i) ~= 0
                start = sum(num_sys(1:i))+1;
                endi = sum(num_sys(1:i)) + num_sys(i+1);
                H_isb(start:endi, sum(isb_status(1:i))) = 1;
            end
        end
    end
end

end