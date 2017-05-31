function [found] = brute_force(yn, y)
% "brute-force" algorithm for sorting multiple spike types (test EVERY 
% possible spike type combination)
%
% INPUTS:
%        yn    = 1 x N signal vector, to be sorted
%        y     = m x N matrix, each row = 1 spike type
%
% OUTPUT: 
%        found = vector of found spike types

[m,n] = size(y); % # of spike types vs # of samples
Mc = 2^m; % total number of possible combos

% create matrix, each row = list of spike types for each combo
C = zeros(Mc,m); % matrix of all possible combos of spike types
k = 0; % counter for # of combinations made
for i = 1:m % find all possible combos of spike types
    add_C = nchoosek(1:m,i); % matrix of possible combos to be added onto end of C
    [a,b] = size(add_C);
    C((k+1):(k+a),1:i) = add_C; % add on next set of possible combos to end of C
    k = k + a;
end

% create matrix, each row = signal vector for each combo
yc = zeros(Mc,n); % new signal matrix, each row = each possible combo of spike type signals
for i = 1:Mc
    for j = 1:m % sum spike shapes in each combo
        if C(i,j) ~= 0, yc(i,:) = yc(i,:) + y(C(i,j),:); end
    end
end

% brute-force algorithm
for i = 1:Mc % test which combo brings signal vector closest to origin
    w(i) = norm(yc(i,:)-yn)^2;
end

[M, I] = min(w);
found = C(I,:); % record spike shape combo that brings signal vector closest to origin
found = found(found>0); % get rid of 0's (null spike shapes)
end