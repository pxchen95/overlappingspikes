function [found] = greedy_sort(yn, y, repl)
% forward greedy algorithm for sorting multiple spike types, with option 
% for with or without replacement
%
% INPUTS:
%        yn    = 1 x N signal vector, to be sorted
%        y     = m x N matrix, each row = 1 spike type
%        repl  = 0 for without replacement, 1 for with replacement
%
% OUTPUT: 
%        found = vector of found spike types

[m,n] = size(y); % # of spike types vs # of samples
found = []; l = 1; % vector (and its counter) for list of detected spike shapes
s = 0; % stop signal for while loop
is_avail = 1:m; % which spike shapes are available to be tested
while s == 0
    for j = is_avail % test which shape brings signal vector closest to origin
        if j ~= 0, w(j) = norm(y(j,:)-yn)^2; end
    end
    w(is_avail == 0) = NaN; 
    [M, I] = min(w);
         
    if M < norm(yn)^2 
        found(l) = I; % record spike shape that brings signal vector closest to origin
        yn = yn - y(I,:); % update signal vector with detected shape removed
        l = l+1; % update counter (# of detected spikes)
        if repl == 0, is_avail(I) = 0; end % w/o replacement - don't test accepeted shape again
    else s = 1; % if no spike shape works, stop searching
    end 
end

end