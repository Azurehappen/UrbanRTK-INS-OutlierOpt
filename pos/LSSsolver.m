function [pos,clock_bias,res] = LSSsolver(p,xk,H_offset,cpt)
 
y = cpt.corr_range;
num = length(y); % The number of measurement
H = zeros(num,4);
R = zeros(num,1); 
r = zeros(num,1);
off = zeros(num,1);
if p.post_mode == 1
    s_pos_ecef = cpt.s_pos_prc;
else
    s_pos_ecef = cpt.s_pos_ecef;
end
Es = zeros(num,1);
for iter=1:p.Nls
    for j=1:num
    R(j)=norm(s_pos_ecef(:,j)-xk(1:3));
    V= (xk(1:3)-s_pos_ecef(:,j))'/R(j)+...
        [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0]; 
    H(j,:)=[V 1]; 
    if p.post_mode == 1
%         Es(j) = norm(cpt.s_pos_prc(:,j)-xk(1:3))+sagnac(p,cpt.s_pos_prc(:,j),xk)...
%             -norm(cpt.s_pos_ecef(:,j)-xk(1:3))-sagnac(p,cpt.s_pos_ecef(:,j),xk);
         Es(j)=0;
    end
    r(j) = R(j)+sagnac(p,s_pos_ecef(:,j),xk);
    if ~isempty(H_offset)
         ind = find(H_offset(j,:)==1);
         if ~isempty(ind)
              off(j) = xk(4+ind);
         end
    end
    end
res = y-Es - r - xk(4)-off;
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