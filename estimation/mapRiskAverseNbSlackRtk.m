function [flag,x_post,P_post,b,J_out,augcost,num_nodes,constraint,risk,penalty_sum] =...
    mapRiskAverseNbSlackRtk(p, num_constrain,y,H,P,R,J_l,x_prior)
% Soft-constrained DiagRAPS for RTK Float solutions.
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

num_states = p.modeToNumUserErrStates(p.state_mode);
[m,n] = size(H);
ind_n = [1:m/3,m*2/3+1:m]; % index without phase
P_n = P(1:num_states,1:num_states);
R_n = R(ind_n,ind_n);
Jpminus_n = P_n^-1;         % state prior info. matrix no phase
Jrminus_n = R_n^-1;         % measurement info. matrix no phase
H_n = H(ind_n,1:num_states);
y_n = y(ind_n);
x_prior_n = x_prior(1:num_states);

[m_n,n_n] = size(H_n);
lowerbound = zeros(m_n,1); % lower bound on measurement selection vector b
upperbound = ones(m_n,1);  % upper bound on measurement selection vector b

% see intlinprog function documentation for details
diagR  = diag(R_n);               % diagonal entries of measurement covariance
diagRs = sqrt(diagR);
sH     = diag(diagRs) \ H_n;      % scale row i by {\sigma_i)^{-1}: diagRs^-1 * H
ieqLHS = -(sH.*sH)';            % inequality constraint LHS matrix (-G in the paper)
ieqLHS = ieqLHS(1:num_constrain,:);        % remove unconstrained rows
ieqRHS = diag(Jpminus_n(1:num_constrain,1:num_constrain)...
    - J_l(1:num_constrain,1:num_constrain));   % inequality constraint RHS vector (-g in the paper)
% ieqRHS = ieqRHS(1:num_constrain,:);         % remove unconstrained rows
option = optimoptions(@linprog,'display','off'); % for output supression
 
constraint = -ieqLHS*ones(m_n,1)+diag(Jpminus_n(1:num_constrain,1:num_constrain));

E_P_n = chol(Jpminus_n);    % cholesky decomp. of state prior info. matrix
E_R_n = chol(Jrminus_n);    % cholesky decomp. of measurement info. matrix

% ensure feasibility
b = ones(m_n,1);          % Initial b to be 1
num_nodes = 0;
flag = true;
ind = [];
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
    lowerbound = [zeros(m_n+length(ind),1)]; % lower bound on measurement selection vector b+s
    upperbound = [ones(m_n,1);upper_s];  % upper bound on measurement selection vector b+s
    flag = false;
end
penalty_sum = 0;
xcurr = x_prior_n;
num_iter = 1;
total_trial = 15;
b = ones(m_n,1);
res    = y_n-H_n*xcurr;
s_res  = res./diagRs;
cost_b = (s_res).*(s_res);
if flag == false
    penalty = 50*ones(length(ind),1);
    cost_b = [cost_b;penalty];
    C(1) = cost(xcurr,x_prior_n,Jpminus_n,[b;upper_s],cost_b);
else
    C(1) = cost(xcurr,x_prior_n,Jpminus_n,b,cost_b);
end

Aeq = [];
Beq = [];
while num_iter < total_trial
    res    = y_n-H_n*xcurr;             % residual
    s_res  = res./diagRs;           % residual scaled by meas. std (z-value)
    cost_b = (s_res).*(s_res);      % cost function coefficients for b optimization
    if flag == false
        penalty = 50*ones(length(ind),1);
        cost_b = [cost_b;penalty];
    end
    % solve for optimal integer measurement selection vector
    % solves the optimization problem in eqn (22) in [1].
    [b_now,~,exitflag,~] = linprog(cost_b,ieqLHS,ieqRHS,Aeq,Beq,lowerbound,upperbound,option);
    b = b_now(1:m_n);
    if length(b_now) ~= m_n
        penalty_sum = 50*sum(b_now(m_n+1:end));
    end
    if (exitflag >= 1)
        % feasible solution is found
        C(2*num_iter) = cost(xcurr,x_prior_n,Jpminus_n,b_now,cost_b); % C(x_{l-1}, b_l)
        [x_post, ~]  = MAP(b,y_n,H_n,E_R_n,E_P_n,x_prior_n); % compute x_l
    else
        % feasible solution not found 
        % Check the feasiblity of G*b before try.
        warning('feasible solution not found, this should not happen.')
    end
    
    xcurr = x_post;
    C(2*num_iter+1) = cost(xcurr,x_prior_n,Jpminus_n,b_now,cost_b); % C(x_l, b_l)
    augcost = C(2*num_iter+1);
    if abs(C(2*num_iter+1)-C(2*num_iter-1))<0.001
        break
    end
    num_iter = num_iter + 1;
end % while cnt_it

Jpminus = P^-1;         % state prior info. matrix
Jrminus = R^-1;         % measurement info. matrix
E_P = chol(Jpminus);    % cholesky decomp. of state prior info. matrix
E_R = chol(Jrminus);    % cholesky decomp. of measurement info. matrix
b_all = [b(1:m_n/2);b(1:m_n/2);b(m_n/2+1:end)];

[x_post, augcost]  = MAP(b_all,y,H,E_R,E_P,x_prior); % compute x_l
risk = augcost;

J_out = calcJb(b_all,H,R,Jpminus); % posterior state information matrix
%J_out % should meet spec
P_post = inv(J_out);

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