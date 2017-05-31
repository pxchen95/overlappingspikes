function [found] = greedy_pairs(yn, y, m, gamma, eta, repl)
% greedy algorithm for sorting multiple spike types with time shifts
% (test every shift and every pair)
%
% INPUTS:
%        yn    = 1 x N signal vector, to be sorted
%        y     = m*ns x N matrix, each row = 1 spike type (incl. shifts as types)
%        m     = number of spike types
%        gamma = probability of each spike type (incl. shifts as types)
%        eta   = noise level, 1/SNR
%
% OUTPUT: 
%        found = vector of found spike types

[mm,nn] = size(y); % # of spike types (incl. shifts as types) vs # of samples
ns = mm/m; % number of time shifts

if m > 1, U = 2; else U = 1; end

Mc = 0; Mt = 0; % initialize
for i = 1:U
    Mc = Mc + nchoosek(m,i)*ns^i; % total number of possible combos (incl. time shifts)
    Mt = Mt + nchoosek(m,i); % total number of possible combos (spike TYPES only)
end

k = 0; % counter for # of combinations made
Ct = zeros(Mt,U); % matrix of all possible combos of spike types
for i = 1:U % find all possible combos of spike types
    add_C = nchoosek(1:m,i); % matrix of possible combos to be added onto end of C
    [a,b] = size(add_C);
    Ct((k+1):(k+a),1:i) = add_C; % add on next set of possible combos to end of C
    k = k + a;
end

% create matrix, each row = list of spike types (incl. shifts as types) for each combo
C = zeros(Mc,U); % matrix of all possible combos of spike types (incl. shifts as types)
for i = 1:m
    types(i,:) = ((i-1)*ns+1):i*ns; % each row = index of time shifts for each spike type
end

% fill in C, matrix of all possible combos of spike types (incl. shifts as types)
l = 1; % counter along C
for i = 1:Mt
    vectors = {}; k = 1; % reset
    for j = 1:U
        if Ct(i,j) ~= 0, vectors{k} = types(Ct(i,j),:); k = k + 1; end
    end
    if ~isempty(vectors)
        combs = combine(vectors); [a,b] = size(combs); 
        C(l:(l+a-1),1:b) = combs; l = l + a; 
    end
end

% create matrix, each row = signal vector for each combo
yc = zeros(Mc,nn); % new signal matrix, each row = each possible combo of spike type signals
for i = 1:Mc
    for j = 1:U
        if C(i,j) ~= 0, 
            yc(i,:) = yc(i,:) + y(C(i,j),:); % sum spike shapes in each combo 
        end
    end
end

% for penalty (adjusted prior probabilities)
prior0 = -2*eta^2*log(1 - gamma); % without spike
prior1 = -2*eta^2*log(gamma/ns); % with spike

% calculate penalties for pairs (gamma)
prior0_C = zeros(Mc, 1); prior1_C = zeros(Mc, 1);
prior0_C(1:mm) = prior0; prior1_C(1:mm) = prior1; 
for i = (mm+1):Mc
    prior0_C(i) = sum(prior0(C(i,:)));
    prior1_C(i) = sum(prior1(C(i,:)));
end
% prior0_C = zeros(Mc, 1); prior1_C = zeros(Mc, 1);

% ---------------- greedy-algorithm ----------------
epsilon = 0; % 0 = accept spikes that cause any improvement
is_avail = 1:Mc; % which spike shapes are available to be tested

found = []; % list of detected spike shapes
k = 1; % # of forward steps (or # of acceptances) + 1
s = 0; % stop signal for FoBa algorithm
lambda = zeros(Mc, 1); % penalty for adding a spike
ynsq(1) = norm(yn)^2; % original error

while s == 0    
    w = NaN*zeros(1,Mc); % initialize
    for j = is_avail % test which shape brings signal vector closest to origin
        if j ~= 0
            w(j) = norm(yc(j,:) - yn)^2 + prior1_C(j) - prior0_C(j); % error term incl. penalty
        end
    end
    [M, I] = min(w);
    ynsq(k+1) = M - prior1_C(I) + prior0_C(I);
    delta_f = ynsq(k) - M; % improvement from adding a spike (to found)

    if delta_f > epsilon
        found(k) = I; % record spike shape that brings signal vector closest to origin
        yn = yn - yc(I,:); % update signal vector with detected shape removed
        k = k + 1;
        if repl == 0, % w/o replacement - don't test accepeted type (or its shifts) again
            for j = C(I,:)
               I_type = ceil(j/ns);
               is_avail(sum(ceil(C/ns) == I_type, 2)>0) = 0;
            end
        end 
    else
        s = 1; % if all spike shapes cause increase in error, stop searching
    end 
end

found = C(found,:); % makes matrix
found = found(:); % makes column vector
found = found(found>0); % kills 0s
end