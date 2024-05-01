function [Phi, Qd] = computeDiscreteModel(F, G, W, dt)
% Compute the Phi and Qd increments
% Parameters:                                                              
%   F:          the continuous time state transition matrix                    
%   G:          the continuous time process noise covariance matrix            
%   W:          the process noise covariance matrix 
%   dt:         the time step duration

[m,n] = size(F);
% pre-compute the common matrix Q.
Q = G*W*G';
% Farrell eq.4.113
chi = [ -F        Q ;
    zeros(m)  F'];

% Farrell eq.4.114
gamma = expm(chi*dt);

% Farrell eq.4.115
Phi = gamma((m+1):(2*m),(m+1):(2*m))';

% Farrell eq.4.116
Qd  = Phi*gamma(1:m,(m+1):(2*m));