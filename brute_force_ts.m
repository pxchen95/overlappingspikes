function [found] = brute_force_ts(yn, y, m, gamma, eta)
% "brute-force" algorithm for sorting multiple spike types with time shifts
% (test every shift, every pair, and every triplet)
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

Mc = 0; Mt = 0; % initialize
for i = 0:m
    Mc = Mc + nchoosek(m,i)*ns^i; % total number of possible combos (incl. time shifts)
    Mt = Mt + nchoosek(m,i); % total number of possible combos (spike TYPES only)
end

k = 0; % counter for # of combinations made
Ct = zeros(Mt,m); % matrix of all possible combos of spike types
for i = 1:m % find all possible combos of spike types
    add_C = nchoosek(1:m,i); % matrix of possible combos to be added onto end of C
    [a,b] = size(add_C);
    Ct((k+1):(k+a),1:i) = add_C; % add on next set of possible combos to end of C
    k = k + a;
end

% create matrix, each row = list of spike types (incl. shifts as types) for each combo
C = zeros(Mc,m); % matrix of all possible combos of spike types (incl. shifts as types)
for i = 1:m
    types(i,:) = ((i-1)*ns+1):i*ns; % each row = index of time shifts for each spike type
end

% fill in C, matrix of all possible combos of spike types (incl. shifts as types)
l = 1; % counter along C
for i = 1:Mt
    vectors = {}; k = 1; % reset
    for j = 1:m
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
    for j = 1:m
        if C(i,j) ~= 0, yc(i,:) = yc(i,:) + y(C(i,j),:); end % sum spike shapes in each combo
    end
end

% for penalty (adjusted prior probabilities)
prior0 = -2*eta^2*log(1 - gamma); % without spike
prior1 = -2*eta^2*log(gamma/ns); % with spike

% calculate penalties for pairs (gamma)
prior0_C = zeros(Mc, 1); prior1_C = zeros(Mc, 1);
prior0_C(1:mm) = prior0; prior1_C(1:mm) = prior1; 
for i = (mm+1):Mc
    prior0_C(i) = sum(prior0(find(C(i,:))));
    prior1_C(i) = sum(prior1(find(C(i,:))));
end
% prior0_C = zeros(Mc, 1); prior1_C = zeros(Mc, 1);

% brute-force algorithm
for i = 1:Mc % test which combo brings signal vector closest to origin
    w(i) = norm(yc(i,:) - yn)^2 + prior1_C(i) - prior0_C(i);
end

[M, I] = min(w);
found = C(I,:); % record spike shape combo that brings signal vector closest to origin
found = found(found>0); % get rid of 0's (null spike shapes)
end