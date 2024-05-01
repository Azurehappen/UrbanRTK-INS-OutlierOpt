function [pos,clock_bias,res,cost] = LTSsolver(p,xk,H_offset,s_pos_ecef,y)

num = length(y); % The number of measurement
H   = zeros(num,4);
R   = zeros(num,1);
r   = zeros(num,1);
off = zeros(num,1);


switch p.linear

case 0 % Nonlinear model
    mu_x = xk;
    for iter=1:p.Nls
        for j=1:num
            R(j)=norm(s_pos_ecef(:,j)-xk(1:3));
            V= (xk(1:3)-s_pos_ecef(:,j))'/R(j);        
            H(j,:)=[V 1];  
            r(j) = R(j)+sagnac(p,s_pos_ecef(:,j),xk);
            if ~isempty(H_offset)
                ind = find(H_offset(j,:)==1);
                off(j) = xk(4+ind);
            end
        end
        res   = y- r - xk(4)-off;
        H_os  = [H,H_offset];
        Hdim2 = size(H_os,2);
        priorposCov = diag(priorposstd(1:Hdim2).^2);
        yCov = diag(p.sig_y^2.*ones(num,1));
        E_R  = sqrt(yCov^(-1));
        E_P  = sqrt(priorposCov^(-1));    
        A    = [E_R*H_os; E_P];
        c    = [E_R*res; 0];
        g    = []; % # of residuals used for Cost minimization
        if isempty(g)
            [rew,~] = ltsregres(A,c,'plots',0,'intercept',0);
        else
            [rew,~] = ltsregres(A,c,'plots',0,'intercept',0,'h',g);
        end
        delta_x = rew.slope;
        b   = rew.flag;
        xk  = xk+delta_x; 
        
        Phi = @(b) diag(b);
        b_y = b(1:num); b_x = b(num+1:end);
        A_b = [E_R*Phi(b_y)*H_os; E_P*Phi(b_x)];
        c_b = [E_R*Phi(b_y)*res; E_P*Phi(b_x)*mu_x];
        
        if nargout == 5
            cost = norm(A_b*xk - c_b); 
        end
        
        if (norm(delta_x) < p.LSthrsh)
            break;
        end
        if (iter>p.Nls)&& (norm(delta_x) > p.LSthrsh)
        warning('Postion path length iteration failed in user_pos calculation');
        end    
    end
    
    

case 1  % Linear model
    
    for j=1:num
        R(j)=norm(s_pos_ecef(:,j)-xk(1:3));
        V = (xk(1:3)-s_pos_ecef(:,j))'/R(j);        
        H(j,:)=[V 1];  
        r(j) = R(j)+sagnac(p,s_pos_ecef(:,j),xk);
        if ~isempty(H_offset)
            ind = find(H_offset(j,:)==1);
            if ~isempty(ind)
                off(j) = xk(4+ind);
            end
        end
    end
    res   = y - r - xk(4)-off;
    H_os  = [H,H_offset];
    Hdim2 = size(H_os,2);
    priorposCov = diag((p.priorposstd(1:Hdim2)).^2);
    yCov  = diag(p.sig_y^2.*ones(num,1));
    E_R   = sqrt(yCov^(-1));
    E_P   = sqrt(priorposCov^(-1));
    mu_x  = xk;
    A     = [E_R*H_os; E_P];
    c     = [E_R*res; E_P*mu_x];
    
    [rew,~] = LTS(A,c);        
    delta_x = rew.slope; 
    b  = rew.flag; % Binary Decision Vector
    xk = xk + delta_x;
    
    Phi = @(b) diag(b);
    b_y = b(1:num); b_x = b(num+1:end);
    A_b = [E_R*Phi(b_y)*H_os; E_P*Phi(b_x)];
    c_b = [E_R*Phi(b_y)*res; E_P*Phi(b_x)*mu_x]; 
    
    if nargout == 5
        cost = norm(A_b*xk - c_b); 
    end

end
pos = xk(1:3);
clock_bias = xk(4);



end