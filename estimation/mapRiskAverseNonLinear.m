function [x_post,P_post,b,augcost,J_out,exitflag, num_iter, comp_t] =...
    mapRiskAverseNonLinear(y,H,P,R,J_l,x_prior,Rot_e2g)
% Solves RAPS using B&B integer optimization approach.
% Computes MAP state estimate using the selected measurements.
% OUTPUT:   x_post   - posterior state vector estimate
%           by       - measurement selection vector (binary)
%           augcost  - augmented cost for RAPS B&B
%           exitflag - see MATLAB function: intlinprog for description
% INPUT:    res - measurements residuals: y - f(x0) = H*dx, where dx =
% x_plus - x_minus
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
ieqLHS = ieqLHS(1:3,:);         % remove unconstrained rows
ieqRHS = diag(Jpminus(1:3,1:3) - J_l);   % inequality constraint RHS vector (-g in the paper)
ieqRHS = ieqRHS(1:3,:);         % remove unconstrained rows
intcon = 1:1:m;                 % entries of b vector that can take only integer values
option = optimoptions(@intlinprog,'display','off'); % for output supression

% ensure feasibility
loosen_count = 0;

E_P = chol(Jpminus);    % cholesky decomp. of state prior info. matrix
E_R = chol(Jrminus);    % cholesky decomp. of measurement info. matrix

dxcurr = zeros(length(x_prior), 1);
num_iter = 1;
total_trial = 10;
b = ones(m,1);
C = zeros(1, 2*total_trial);
C(1) = cost(dxcurr,P,b,H,y,R);
tic
while num_iter < total_trial
    Aeq = [];   % Coefficients for equality constrains: Aeq * b = Beq
    Beq = [];   % Values for equality constrains: Aeq * b = Beq
    res    = y-H*dxcurr;             % residual
    s_res  = res./diagRs;           % residual scaled by meas. std (z-value)
    cost_b = (s_res).*(s_res);      % cost function coefficients for b optimization

    % solve for optimal integer measurement selection vector (Branch & Bound search)
    % solves the optimization problem in eqn (22) in [1].
    oldb = b;
    [b,~,exitflag,~] = intlinprog(cost_b,intcon,ieqLHS,ieqRHS,Aeq,Beq,lowerbound,upperbound,option);
    switch exitflag
        case 3
            disp('Optimal solution found with poor constraint feasibility.')
        case 2
            disp('Solver stopped prematurely. Integer feasible point found.')
        case 1
            %disp('Optimal solution found.')
        case 0
            disp('Solver stopped prematurely. No integer feasible point found.')
        case -1
            disp('Solver stopped by an output function or plot function.')
        case -2
            disp('No feasible point found.')
        case -3
            disp('Root LP problem is unbounded.')
        case -9
            % This case tends to occur when the ratio of the largest to lowest cost is too great.
            % Our solution is to constrain the problem to exclude the
            % highest cost measurement to be excluded, until the issue
            % goes away.
            [srt_cst,ind_srt] = sort(cost_b);
            neqcnstrt = 0;  % number of equality constraints
            while exitflag == -9 && neqcnstrt < m/2
                disp('Solver lost feasibility probably due to ill-conditioned cost. Removing the largest cost')
                index_max = find(cost_b == srt_cst(m-neqcnstrt));
                cnt = length(index_max);
                for i=1:cnt
                    Aeq(neqcnstrt+1,:) = zeros(1,m);  % initialize new row
                    Aeq(neqcnstrt+1,index_max(i)) = 1;% constrain the correct element
                    Beq(neqcnstrt+1,1) = 0;           % constrain to be excluded
                end
                neqcnstrt = neqcnstrt + cnt;
                [b,~,exitflag,~] = intlinprog(cost_b,intcon,ieqLHS,ieqRHS,Aeq,Beq,lowerbound,upperbound,option);
                if (exitflag <= 0)
                    disp('Solver still unsuccessful')
                    exitflag
                end
            end
        otherwise
            exitflag
            disp('unknown case')
    end

    if (exitflag >= 1)
        % feasible solution is found
        C(2*num_iter) = cost(dxcurr,P,b,H,y,R); % C(x_{l-1}, b_l)
        [x_post, augcost, dx_post]  = MAP(b,y,H,E_R,E_P,x_prior,Rot_e2g); % compute x_l
        flg_feasbl = true;
    else
        % feasible solution not found
        fprintf('\nExitflag = %1.0f. Feasible solution not found.\n',exitflag)
        flg_feasbl = false;
        if (loosen_count > 3)
            break;
        end
        b = ones(m,1);
        max_neg_lhs = ieqLHS * b; 
        ind = find(ieqRHS < max_neg_lhs);
        ieqRHS(ind) = 0.8 * max_neg_lhs(ind);    % loosen constraint
        dx_post = dxcurr;                          % try again
        loosen_count = loosen_count + 1;
        continue;
    end
    %b'
    dxcurr = dx_post;
    C(2*num_iter+1) = cost(dxcurr,P,b,H,y,R); % C(x_l, b_l)
    if and(abs(C(2*num_iter+1)-C(2*num_iter-1))<0.001,flg_feasbl)
        break
    end
    num_iter = num_iter + 1;
end % while cnt_it
if ~flg_feasbl
    % if still not found feasible solution, no RAPS performs.
    disp('Feasible solution not found, perform general MAP.')
    b = ones(m,1);
    [x_post, augcost,~]  = MAP(b,y,H,E_R,E_P,x_prior, Rot_e2g);
end
J_out   = calcJb(b,H,R,Jpminus); % posterior state information matrix
%J_out % should meet spec
P_post = inv(J_out);

comp_t = toc; % computation time

end

function [C] = cost(dx,P,b,H,y,R)
Pi = inv(P);
ey = y - H*dx;
Ri = inv(R);
Pb = diag(b);
if isempty(b)
    C1 = dx' * Pi * dx;
    C2 =  ey' *  Ri *  ey;
    C  = C1 + C2;
else
    C = dx' * Pi * dx + ey' * Pb * Ri * Pb * ey;
end
end

%--------------------------------------------------------------------------
function [x_post,aug_cost,dx] = MAP(by,y,H,E_R,E_P,x_prior, Rot_e2g)
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
    c = [E_R * Phiby * y; zeros(length(x_prior),1)];   % eqn. (10) in [1]
    dx = (A'*A)^-1*A'*c; % LeastSquares % posterior state estimate
    aug_cost = norm(A*dx-c)^2;          % eqn. (11) in [1]
end
x_post = x_prior + Rot_e2g' * dx;

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