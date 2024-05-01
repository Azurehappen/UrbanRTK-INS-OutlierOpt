function [pos,clock_bias,res] = userpos_2diff(p,cpt)
% This is solver for computing user position, satellite system offset.
% Using double difference. Receiver clock terms been taken off
% Measurement selection applied
% Input: 
%       s_pos_ecef: 3-by-N Satellite position in ECEF frame.
%       x0 : 3-by-1 initial interative coordinates in ECEF frame.
%       y: m-by-1 Corrected pseudorange.
%
% Output:
%       
%       
%       
%       

%-------------------%
% Initialize
y = cpt.corr_range;
num_sys = cpt.num_sv;
if p.post_mode == 1
    if p.IGS_enable == 1
        sat_R = cpt.s_pos_prc;
    else
        sat_R = cpt.s_pos_ecef;
    end
else
    sat_R = cpt.s_pos_ecef;
end
x0 = p.state0(1:3);
[H_offset,x_offset] = sys_offset(num_sys);
xk = [x0;x_offset];
%------------------%
pivot = find(cpt.elev(1:cpt.num_sv(1))==max(cpt.elev(1:cpt.num_sv(1)))); %Find the pivot satellite
pivot_R = y(pivot); y(pivot) = [];
pivot_sat = sat_R(:,pivot); sat_R(:,pivot)=[];
num = length(y); % The number of measurement
H = zeros(num,3);
R = zeros(num,1);
r = zeros(num,1);
off = zeros(num,1);
Es = zeros(num,1);
y = y - pivot_R;
for iter=1:p.Nls
    R0 = norm(pivot_sat - xk(1:3));
    for j=1:num
    R(j)=norm(sat_R(:,j)-xk(1:3));
    H(j,:) = (xk(1:3)-sat_R(:,j))'/R(j)-(xk(1:3)-pivot_sat)'/R0;        
    if p.post_mode == 1
         Es(j) = 0;
%          Es(j) = norm(cpt.sat_posprc_Rcorr(:,j)-xk(1:3))...
%              -norm(cpt.sat_pos_Rcorr(:,j)-xk(1:3));
    end
    r(j)=R(j)+sagnac(p,sat_R(:,j),xk)-R0 - sagnac(p,pivot_sat,xk) ;

    if ~isempty(H_offset)
         ind = find(H_offset(j,:)==1);
         if ~isempty(ind)
              off(j) = xk(4+ind);
         end
    end
    end
res = y - Es - r -off;
H_os = [H,H_offset];
delta_x = (H_os'*H_os)^(-1)*H_os'*(res);
xk=xk+delta_x; 
if (norm(delta_x) < p.LSthrsh)
    break;
end
if (iter>p.Nls)&& (norm(delta_x) > p.LSthrsh)
     warning('Postion path length iteration failed in user_pos calculation');
end  

end
GDOP = sqrt(trace((H_os'*H_os)^(-1)));
if GDOP_check(p,GDOP)==1  
    pos = xk(1:3);
    clock_bias = 0;
else
    pos = [];
    clock_bias = [];
end
res = [res(1:pivot-1);0;res(pivot:end)];


end
