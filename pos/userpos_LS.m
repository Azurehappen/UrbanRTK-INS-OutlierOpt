function [pos,clock_bias,res] = userpos_LS(p,s_pos_ecef,x0,y,num_sys)
% This is a least square solver for computing user position, receiver clock
% bias, satellite system offset.
%
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
num = length(y); % The number of measurement
[H_offset,x_offset] = sys_offset(num_sys);
H = zeros(num,4);
R = zeros(num,1);
r = zeros(num,1);
off = zeros(num,1);
xk = [x0;x_offset];
%------------------%
for iter=1:p.Nls
    for j=1:num
        R(j)=norm(s_pos_ecef(:,j)-xk(1:3));
        V= (xk(1:3)-s_pos_ecef(:,j))'/R(j);        
        H(j,:)=[V 1];  
        r(j) = R(j)+sagnac(p,s_pos_ecef(:,j),xk);
        if ~isempty(H_offset)
            ind = find(H_offset(j,:)==1);
            if ~isempty(ind)
                off(j) = xk(4+ind);
            end
        end
    end
    res = y- r - xk(4)-off;
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
pos = xk(1:3);
clock_bias = xk(4);
end
