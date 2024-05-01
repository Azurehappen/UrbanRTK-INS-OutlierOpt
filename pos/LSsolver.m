function [pos,clock_bias,isb,res] = LSsolver(p,xk,H_isb,cpt)
 
y = cpt.corr_range;
num = length(y); % The number of measurement
H = zeros(num,4);
R = zeros(num,1); 
r = zeros(num,1);
off = zeros(num,1);
if p.post_mode == 1
    if p.IGS_enable == 1
        s_pos_ecef = cpt.s_pos_prc;
    else
        s_pos_ecef = cpt.s_pos_ecef;
    end
else
    s_pos_ecef = cpt.s_pos_ecef;
end
% if p.mk == 1
%     xk = [xk;0]; %Add Iono factor term for PPP
% end
for iter=1:p.Nls
    for j=1:num
    R(j)=norm(s_pos_ecef(:,j)-xk(1:3));
    V= (xk(1:3)-s_pos_ecef(:,j))'/R(j)+...
        [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0]; 
    H(j,:)=[V 1];
    r(j) = R(j)+sagnac(p,s_pos_ecef(:,j),xk);
    if ~isempty(H_isb)
         ind = find(H_isb(j,:)==1);
         if ~isempty(ind)
              off(j) = xk(4+ind);
         end
    end
    end
%     if p.mk==1
%         res = y - r - xk(4)-off - cpt.IoFac*xk(end);
%         H_os = [H,H_offset,cpt.IoFac];
%     else
        res = y - r - xk(4)-off;
        H_os = [H,H_isb];
%     end
delta_x = (H_os'*H_os)^(-1)*H_os'*(res);
xk=xk+delta_x; 
if (norm(delta_x) < p.LSthrsh)
     break;
end 
if (iter>p.Nls)&& (norm(delta_x) > p.LSthrsh)
     warning('Postion path length iteration failed in user_pos calculation');
end  

end
% pos = xk(1:3);
% clock_bias = xk(4);
GDOP = sqrt(trace((H_os'*H_os)^(-1)));
isb = [];
if GDOP_check(p,GDOP)==1  
    pos = xk(1:3);
    clock_bias = xk(4);
    if ~isempty(H_isb)
        isb = xk(5:end);
    end
else
    pos = [];
    clock_bias = [];
end      

end