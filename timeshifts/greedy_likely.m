function [found] = greedy_likely(yn, y, ns, eta, gamma, repl)
% forward greedy algorithm WITH PENALTY for detecting multiple spikes,
% option for with or without replacement
%
% INPUTS:
%        yn    = 1 x N signal vector, to be sorted
%        y     = m x N matrix, each row = 1 spike type (incl. shifts as types)
%        ns    = number of time shifts
%        eta   = noise level, 1/SNR
%        gamma = probability of each spike type (incl. shifts as types)
%        repl  = 0 for without replacement, 1 for with replacement
%
% OUTPUT: 
%        found = vector of found spike types

epsilon = 0; % 0 = accept spikes that cause any improvement

[m,n] = size(y); % # of spike types vs # of samples
is_avail = 1:m; % which spike shapes are available to be tested

found = []; % list of detected spike shapes
k = 1; % # of forward steps (or # detected spikes) + 1

% for penalty (adjusted prior probabilities)
prior0 = -2*eta^2*log(1 - gamma); % without spike
prior1 = -2*eta^2*log(gamma/ns); % with spike

ynsq(1) = norm(yn)^2; % original error
s = 0; % stop signal
while s == 0
    w = NaN*zeros(1,m);
    for j = is_avail % test which shape brings signal vector closest to origin
        if j ~= 0
            w(j) = norm(y(j,:) - yn)^2 + prior1(j) - prior0(j); % error term
        end
    end
    [M, I] = min(w);
    ynsq(k+1) = M - prior1(I) + prior0(I);
    delta_f = ynsq(k) - M; % improvement from adding a spike (to found)
    
    if delta_f > epsilon % accept proposed spike type
        found(k) = I; % record spike shape that brings signal vector closest to origin
        yn = yn - y(I,:); % update signal vector with detected shape removed
        k = k + 1;
        if repl == 0, % w/o replacement - don't test accepeted type (or its shifts) again
            I_type = (1:ns) + (ceil(I/ns)-1)*ns; 
            is_avail(I_type) = 0; 
        end 
    else
        s = 1; % if all spike shapes cause increase in error, stop searching
    end 
end

end