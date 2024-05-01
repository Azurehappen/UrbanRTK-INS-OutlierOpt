function [xEstimate,Weights] = Mestimator(y,H,bisqk)

    iterTol   = 1e-4;   % IRWLS iteration tolerance
    TotalIter = 0;      % total iterations before convergence
    BisqK     = 4.685;  % Bisquare Constant
    factor    = 0.675;  % MAD normalization factor
    
    if nargin == 3
        BisqK = bisqk;
    end    
    
    % Caclulate Ordinary LS
    calcOLS   = @(y,H) (H'*H)\H' * y;
    xEst_old  = calcOLS(y,H); % inital estimate of x parameter    
    diff      = inf;          % norm(xestimate_k - xestimates_k-1)
    iterCount = 0;
    calcSigEst = @(res,factor) median(abs(res - median(res))) / factor;
    
    %Iteratively Reweighted Least Squares method 
    while diff >= iterTol & iterCount < 1000       
        res    = y - H*xEst;
        sigEst = calcSigEst(res,factor);
        u      = res./sigEst;       
        m      = size(H,1);
        w      = zeros(m,1);        
        for i = 1:1:m
            %Caluclate weights 
            if abs(u(i)) <= BisqK
                w(i) = (1 - (u(i)/BisqK)^2)^2;
            else
                w(i) = 0;
            end
        end
        W    = diag(w);                  % matrix of weights     
        xEst = inv(H'*W*H)*(H'*W*y);     % Weighted LS estimate of x
        
        diff      = norm(xEst - xEst_old);
        xEst_old  = xEst;
        iterCount = iterCount + 1; 
    end  
    xEstimate = xEst;
    Weights   = W;
    TotalIter = iterCount;
end

