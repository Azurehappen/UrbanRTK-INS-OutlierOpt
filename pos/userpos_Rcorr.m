function [pos,clock_bias,res] = userpos_Rcorr(p,cpt)
% This is solver for computing user position, receiver clock
% bias, satellite system offset.
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
        sat_R = cpt.sat_posprc_Rcorr;
    else
        sat_R = cpt.sat_pos_Rcorr;
    end
else
    sat_R = cpt.sat_pos_Rcorr;
end
x0 = p.state0;
[H_offset,x_offset] = sys_offset(num_sys);
xk = [x0;x_offset];
%------------------%
num = length(y); % The number of measurement
H = zeros(num,4);
R = zeros(num,1);
r = zeros(num,1);
off = zeros(num,1);
Es = zeros(num,1);
for iter=1:p.Nls
    for j=1:num
    R(j)=norm(sat_R(:,j)-xk(1:3));
    V= (xk(1:3)-sat_R(:,j))'/R(j);        
    H(j,:)=[V 1];
    if p.post_mode == 1
         Es(j) = 0;
%          Es(j) = norm(cpt.sat_posprc_Rcorr(:,j)-xk(1:3))...
%              -norm(cpt.sat_pos_Rcorr(:,j)-xk(1:3));
    end
    r(j)=R(j);

    if ~isempty(H_offset)
         ind = find(H_offset(j,:)==1);
         if ~isempty(ind)
              off(j) = xk(4+ind);
         end
    end
    end
res = y - Es - r - xk(4)-off;
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
    clock_bias = xk(4);
else
    pos = [];
    clock_bias = [];
end



end
