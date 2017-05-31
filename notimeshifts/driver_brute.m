% tests brute_force.m; no noise (so should show 100% accuracy)

clear

%---------------------spike shapes-----------------------------
N = 20; % total number of sample times
f1 = @(t, w) exp(-t.^2/(2*(w/7)^2)); % 1-peak spike
y(1,:) = synth(10, 1, 1, 1, f1, 0, N); % spike type 1
y(2,:) = synth(10, 1, 3, 1, f1, 0, N); % spike type 2 (wider spike)
y(3,:) = synth(10, -0.5, 1, 1, f1, 0, N); % upside down spike
[m,n] = size(y); % # of spike types vs # of samples
gamma = ones(1,m)*0.5; % probability of each spike type
%--------------------------------------------------------------

for i = 1:1000 % 1000 trials
    % pick actual spikes to overlap
    pick = rand(size(gamma))<gamma;
    P = 1:m;
    act_list{i} = P(pick==1);

    if sum(pick) > 1 % create signal vector
        yn = sum(y(act_list{i},:)); 
    elseif sum(pick) == 0
        yn = 0;
    else
        yn = y(act_list{i},:);
    end
    
    % detect spikes using brute-force algorithm
    found_list{i} = brute_force(yn,y);
end

% for 100% accuracy, should have 1 for tn and C, 0 for everything else
[tn, fn, C, fp, W] = count(act_list, found_list)