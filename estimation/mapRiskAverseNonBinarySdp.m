function [flag,x_post,P_post,b,J_out,augcost,num_nodes,constraint,risk,penalty_sum] =...
    mapRiskAverseNonBinarySdp(num_constrain,y,H,P,R,J_l,x_prior)
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

num_constrain = 9;

diagR  = diag(R);               % diagonal entries of measurement covariance
diagRs = sqrt(diagR);
sH     = diag(diagRs) \ H;      % scale row i by {\sigma_i)^{-1}: diagRs^-1 * H
ieqRHS = J_l - Jpminus;
ieqRHS = ieqRHS(1:num_constrain,1:num_constrain);         % remove unconstrained rows
J_l_diag = diag(J_l(1:num_constrain,1:num_constrain));
Jpmin_reduce = Jpminus(1:num_constrain,1:num_constrain);
ieqLHS_b = [];
for i=1:m
    temp = sH(i,:)'*sH(i,:);
    ieqLHS_b(:,:,i) = temp(1:num_constrain,1:num_constrain);
end

% constraint = -ieqLHS*ones(m,1)+diag(Jpminus(1:num_constrain,1:num_constrain));
constraint = 0;

E_P = chol(Jpminus);    % cholesky decomp. of state prior info. matrix
E_R = chol(Jrminus);    % cholesky decomp. of measurement info. matrix

num_nodes = 0;
penalty_sum = 0;
flag = true;
use_slack = ~isPositiveDefinite(sum(ieqLHS_b,3)-ieqRHS);
if use_slack
    flag = false;
end
u_coe = 500*ones(1,num_constrain);
xcurr = x_prior;
num_iter = 1;
total_trial = 15;
b = ones(m,1); % Initial b to be 1
C = zeros(1, 2*total_trial);
res    = y-H*xcurr;             % residual
s_res  = res./diagRs;           % residual scaled by meas. std (z-value)
cost_b = (s_res).*(s_res);
C(1) = cost(xcurr,x_prior,Jpminus,b,cost_b)+u_coe*J_l_diag;


while num_iter < total_trial
    res    = y-H*xcurr;             % residual
    s_res  = res./diagRs;           % residual scaled by meas. std (z-value)
    cost_b = (s_res).*(s_res);      % cost function coefficients for b optimization
    
    if use_slack
        cvx_begin sdp
            variable b_l(m)
            variable u(num_constrain)
            minimize(cost_b'*b_l + u_coe*u)
            subject to
                sum_Fb = zeros(num_constrain, num_constrain);
                for i = 1:m
                    sum_Fb = sum_Fb + ieqLHS_b(:,:,i) * b_l(i);
                end
                sum_Fb + Jpmin_reduce - diag(J_l_diag-u) == semidefinite(num_constrain);
                %sum_Fb + Jpminus >= diag(J_l_diag-u)
                0 <= b_l <= 1
                0 <= u <= J_l_diag
                % for i = 1:6
                %     0 <= u(i) <= J_l_diag
                % end
                % for i = 7:n
                %     u(i) = 0
                % end
        cvx_end
        penalty_sum = u_coe*u;
    else
        cvx_begin quiet
            variable b_l(m)
            minimize(cost_b'*b_l)
            subject to
                sum_Fb = zeros(num_constrain, num_constrain);
                for i = 1:m
                    sum_Fb = sum_Fb + ieqLHS_b(:,:,i) * b_l(i);
                end
                sum_Fb - ieqRHS == semidefinite(num_constrain);
                0 <= b_l <= 1
        cvx_end
    end

    b=b_l;
    C(2*num_iter) = cost(xcurr,x_prior,Jpminus,b,cost_b)+penalty_sum; % C(x_{l-1}, b_l)
    [x_post, augcost]  = MAP(b,y,H,E_R,E_P,x_prior); % compute x_l
    
    xcurr = x_post;
    C(2*num_iter+1) = cost(xcurr,x_prior,Jpminus,b,cost_b)+penalty_sum; % C(x_l, b_l)
    if abs(C(2*num_iter+1)-C(2*num_iter-1))<0.001
        break
    end
    num_iter = num_iter + 1;
end % while cnt_it
risk = compute_risk(xcurr,x_prior,Jpminus,b,H,y,Jrminus);
J_out   = calcJb(b,H,R,Jpminus); % posterior state information matrix
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

function isPSD = isPositiveDefinite(A)
    % Check if the matrix is square
    if ~isequal(size(A,1), size(A,2))
        error('Matrix must be square.');
    end

    % Compute eigenvalues
    eigenvalues = eig(A);
    ind = eigenvalues<0 & eigenvalues >=-0.0001;
    eigenvalues(ind) = 0;

    % Check if all eigenvalues are non-negative
    isPSD = all(eigenvalues >= 0); % Using a small threshold to account for numerical errors
end