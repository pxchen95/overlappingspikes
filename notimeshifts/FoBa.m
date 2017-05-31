function [found] = FoBa(yn, y, repl)
% Forward-Backward greedy algorithm 
%
% INPUTS:
%        yn    = 1 x N signal vector, to be sorted
%        y     = m x N matrix, each row = 1 spike type
%        repl  = 0 for without replacement, 1 for with replacement
%
% OUTPUT: 
%        found = vector of found spike types

epsilon = 0; % 0 = accept spikes that cause any improvement
nu = 0.5; % parameter for how much worse a backward step can get in order to be accepted

[m,n] = size(y); % # of spike types vs # of samples
is_avail = 1:m; % which spike shapes are available to be tested

found = []; % list of detected spike shapes
k = 1; % # of forward steps (or # detected spikes) + 1
s = 0; % stop signal for FoBa algorithm
M_old = norm(yn)^2; % original least squares error 

while s == 0
    % forward greedy step
    w = NaN*zeros(1,m);
    for j = is_avail % test which shape brings signal vector closest to origin
        if j ~= 0, w(j) = norm(y(j,:) - yn)^2; end
    end
    [M, I] = min(w);
    delta_f(k) = norm(yn)^2 - M; % improvement from adding a spike (to found)

    if delta_f(k) > epsilon
        found(k) = I; % record spike shape that brings signal vector closest to origin
        yn = yn - y(I,:); % update signal vector with detected shape removed
        k = k + 1;
        if repl == 0, is_avail(I) = 0; end % w/o replacement - don't test accepeted shape again
    else
        break % if all spike shapes cause increase in error, stop searching
    end 

   % backward greedy step
   sb = 0; % stop signal for backward greedy
   while sb == 0 % no backward steps if no forward steps
        W = NaN*zeros(1,m); % don't create zero entries
        for j = found % adding which shape brings yn closest to origin?
            W(j) = norm(yn + y(j,:))^2; 
        end
        [mm, ii] = min(W);
        delta_b = mm - norm(yn)^2; % how much worse the backward step gets from adding spike back 

        if delta_b <= nu*delta_f(k-1)
            yn = yn + y(ii,:); % update yn, put removed shape back
            found(found == ii) = []; % shape added back to yn is removed from detected list
            k = k - 1; % remove one forward step from count
            if repl == 0, is_avail(ii) = ii; end % w/o replacement, can test shape again in future forward steps
        else sb = 1; % if no spike shape works, stop searching
        end 
    end
end

end