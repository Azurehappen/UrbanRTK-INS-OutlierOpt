function [flag,x_post,P_post,b,J_out,augcost,num_iter,constraint,risk,penalty_sum] =...
    mapRiskAverseSlack(num_constrain,y,H,P,R,J_l,x_prior)
% Solves RAPS using B&B integer optimization approach.
% Computes MAP state estimate using the selected measurements.
% OUTPUT:   x_post   - posterior state vector estimate
%           by       - measurement selection vector (binary)
%           augcost  - augmented cost for RAPS B&B
%           exitflag - see MATLAB function: intlinprog for description
% INPUT:    y - measurements
%           H - measurement matrix
%           P - Prior Covariance matrix
%           r - Measurement Covariance matrix
%           J_l - Information Matrix Lower Bound
%           x_prior - prior state vector estimate
% reference: [1] - PPP_RAPS_Linear.pdf

Jpminus = P^-1;         % state prior info. matrix
Jrminus = R^-1;         % measurement info. matrix

[m,n] = size(H);
lowerbound = zeros(m,1); % lower bound on measurement selection vector b
upperbound = ones(m,1);  % upper bound on measurement selection vector b

% see intlinprog function documentation for details
diagR  = diag(R);               % diagonal entries of measurement covariance
diagRs = sqrt(diagR);
sH     = diag(diagRs) \ H;      % scale row i by {\sigma_i)^{-1}: diagRs^-1 * H
ieqLHS = -(sH.*sH)';            % inequality constraint LHS matrix (-G in the paper)
ieqLHS = ieqLHS(1:num_constrain,:);         % remove unconstrained rows
ieqRHS = diag(Jpminus - J_l);   % inequality constraint RHS vector (-g in the paper)
ieqRHS = ieqRHS(1:num_constrain,:);         % remove unconstrained rows
intcon = 1:1:m;                 % entries of b vector that can take only integer values
option = optimoptions(@intlinprog,'display','off'); % for output supression

constraint = -ieqLHS*ones(m,1)+diag(Jpminus(1:num_constrain,1:num_constrain));

E_P = chol(Jpminus);    % cholesky decomp. of state prior info. matrix
E_R = chol(Jrminus);    % cholesky decomp. of measurement info. matrix

% ensure feasibility
b = ones(m,1);          % Initial b to be 1
num_nodes = 0;
flag = true;
if ~all(ieqLHS*b <= ieqRHS)
    % Add slack variable for soft constraint
    max_neg_lhs = ieqLHS * b;
    ind = find(ieqRHS < max_neg_lhs);
    slack_mat = zeros(num_constrain,length(ind));
    for i=1:length(ind)
        % G*b + Jp >= (G b_max + Jp) - s
        % -G*b - s <= -G b_max
        slack_mat(ind(i),i) = -1; % -G*b - s <= -G b_max
    end
    upper_s = -max_neg_lhs(ind);
    ieqRHS(ind) = max_neg_lhs(ind);
    % 0 <= s <= G b_max
    ieqLHS = [ieqLHS, slack_mat];
    lowerbound = [zeros(m+length(ind),1)]; % lower bound on measurement selection vector b+s
    upperbound = [ones(m,1);upper_s];  % upper bound on measurement selection vector b+s
    flag = false;
end
penalty_sum = 0;
xcurr = x_prior;
num_iter = 1;
total_trial = 15;
b = ones(m,1);
C = zeros(1, 2*total_trial);
res    = y-H*xcurr;
s_res  = res./diagRs;
cost_b = (s_res).*(s_res);
if flag == false
    penalty = 50*ones(length(ind),1);
    cost_b = [cost_b;penalty];
    C(1) = cost(xcurr,x_prior,Jpminus,[b;upper_s],cost_b);
else
    C(1) = cost(xcurr,x_prior,Jpminus,b,cost_b);
end

while num_iter < total_trial
    Aeq = [];   % Coefficients for equality constrains: Aeq * b = Beq
    Beq = [];   % Values for equality constrains: Aeq * b = Beq
    res    = y-H*xcurr;             % residual
    s_res  = res./diagRs;           % residual scaled by meas. std (z-value)
    cost_b = (s_res).*(s_res);      % cost function coefficients for b optimization
    if flag == false
        % penalty = 10*ones(length(ieqRHS),1);
        penalty = 50*ones(length(ind),1);
        cost_b = [cost_b;penalty];
    end
    % solve for optimal integer measurement selection vector (Branch & Bound search)
    % solves the optimization problem in eqn (22) in [1].
    [b_all,~,exitflag,output] = intlinprog(cost_b,intcon,ieqLHS,ieqRHS,Aeq,Beq,lowerbound,upperbound,option);
    b = b_all(1:m);
    if length(b_all) ~= m
        penalty_sum = 50*sum(b_all(m+1:end));
    end
    if (exitflag >= 1)
        % feasible solution is found
        num_nodes = num_nodes + output.numnodes;
        C(2*num_iter) = cost(xcurr,x_prior,Jpminus,b_all,cost_b); % C(x_{l-1}, b_l)
        [x_post, ~]  = MAP(b,y,H,E_R,E_P,x_prior); % compute x_l
    else
        % feasible solution not found 
        % Check the feasiblity of G*b before try.
        warning('feasible solution not found, this should not happen.')
    end
    
    xcurr = x_post;
    C(2*num_iter+1) = cost(xcurr,x_prior,Jpminus,b_all,cost_b); % C(x_l, b_l)
    augcost = C(2*num_iter+1);
    if abs(C(2*num_iter+1)-C(2*num_iter-1))<0.001
        break
    end
    num_iter = num_iter + 1;
end % while cnt_it

risk = compute_risk(xcurr,x_prior,Jpminus,b,H,y,Jrminus);

J_out   = calcJb(b,H,R,Jpminus); % posterior state information matrix
%J_out % should meet spec
P_post = inv(J_out);

% flag_pos = true;
% flag_vel = true;
% if flag == true
%     [flag_pos,flag_vel] = rapsValidation(td_lambda+1,P,y,H,R,b,num_range);
% end
% if flag_pos == false || flag_vel == false
%     num_iter = 0;
%     b = thresholdTest(td_lambda,P, y, H, R);
%     [x_post,augcost] = MAP(b,y,H,E_R,E_P,x_prior);
%     J_out   = calcJb(b,H,R,Jpminus); % posterior state information matrix
%     P_post = inv(J_out);
%     flag = false;
% end

end

function [C] = cost(x,prior,Jpminus,b,coe_b)
ex = x-prior;
C = ex' * Jpminus * ex + coe_b'*b;
end

function [C] = compute_risk(x,prior,Jpminus,b,H,y,Jrminus)
ex = x-prior;
ey = y - H*x;
Pb = diag(b);
C = ex' * Jpminus * ex + ey' * Pb * Jrminus * Pb * ey;
end

%--------------------------------------------------------------------------
function [x_post,aug_cost] = MAP(by,y,H,E_R,E_P,x_prior)
% Maximum A Posteriori state estimate and the augmented cost function
% when doing measurement selection
% INPUT: by  - measurement selection vector
%        E_P - sqrt(P^-1)
%        E_R - sqrt(R^-1)

if isempty(by)
    error('Input arg "by" is empty.');
else
    Phiby = diag(by);
    A = [E_R * Phiby * H; E_P];             % eqn. (10) in [1]
    c = [E_R * Phiby * y; E_P * x_prior];   % eqn. (10) in [1]
    x_post = (A'*A)^-1*A'*c; % LeastSquares % posterior state estimate
    aug_cost = norm(A*x_post-c)^2;          % eqn. (11) in [1]
end
end
%--------------------------------------------------------------------------
function J = calcJb(by,H,R,Jpminus)
% Calculates posterior state information matrix when doing measurement selection
% INPUT: by - measurement selection vector
%        Jpminus - State information matrix prior, Jpminus = Pminus^-1

Pby  = diag(by); % Phi(by)
PhiH = Pby * H;
J    = PhiH' * R^-1 * PhiH + Jpminus; % eqn. (13) in [1]
end