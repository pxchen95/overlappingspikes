% runs algorithms with and without time shifts

clear

job = 2; % 1 = without time shifts, 2 = with time shifts
penalty = 0; % 1 = with penalty, 0 = without penalty
noise = 1; repl = 0;
rng(0); % fix seed

% --------------- applies for all jobs --------------- %

x = 0:0.1:0.5; % range of noise levels tested
runs = 200;    %
numii = 5;     %  # of runs per noise level = runs * numii
num_eta = length(x); % number of noise levels tested
k = 1; % counter for number of noise levels tested
l = 1; % counter for waitbar
j = 1; % counter for total number of runs
h = waitbar(0,'Please wait...'); % make wait bar
SAD = []; SAD2 = []; SAD3 = []; % initialize, standard absolute deviation

% ---------------------- begin cases ----------------- %

if job == 1
    %--------------------------------------------------------------
    N = 20; % total number of sample times
    f1 = @(t,w) exp(-t.^2/(2*(w/7)^2)); % 1-peak spike
    f2 = @(t,w) 6*t.*exp(-.5*t.^2/(2*(w/7)^2)); % 2-peak spike
    y(1,:) = synth(10, 1, 1, 1, f1, 0, N); % spike type 1
    y(2,:) = synth(10, 1, 3, 1, f1, 0, N); % spike type 2 (wider spike)
    y(3,:) = synth(10, -.5, 1, 1, f1, 0, N); % upside down, short spike
    y(4,:) = synth(10, 1, 1, 1, f2, 0, N); % down-up spike
    [m,n] = size(y); % # of spike types vs # of samples
    gamma = [0.5 0.3 0.1 0.05];
    %--------------------------------------------------------------
    for eta = x % signal to noise ratio
        for ii = 1:numii
            for i = 1:runs % runs per noise value
                % pick spike shapes present
                pick = rand(size(gamma))<gamma;
                P = 1:m;
                act_list{i} = P(pick==1); % cell array, each element = list of actual spike shapes

                % add normally distributed random noise, proportional to eta
                if noise == 0, eta = 0; end
                yn = sum(y(act_list{i},:),1) + eta*randn(1,N);

                % call algorithms
                if penalty == 0
                    found_list{i} = greedy_pairs(yn, y, m, gamma, 0, repl); % cell array, each element = list of detected spike shapes
                    found_list2{i} = greedy_sort(yn, y, repl); % same as greedy_likely with eta = 0
                    found_list3{i} = brute_force(yn, y); % same as brute_force_ts with eta = 0
                elseif penalty == 1
                    found_list{i} = greedy_pairs(yn, y, m, gamma, eta, repl); % cell array, each element = list of detected spike shapes
                    found_list2{i} = greedy_likely(yn, y, 1, eta, gamma, repl);
                    found_list3{i} = brute_force_ts(yn, y, m, gamma, eta);
                end

                l = l + 1; % update counter for waitbar
                waitbar(l/(num_eta*runs*numii)) % update waitbar
            end
            [tn{k}(ii), fn{k}(ii), C{k}(ii), fp{k}(ii), W{k}(ii)] = count(act_list, found_list, 0); % calculate error fractions
            [tn2{k}(ii), fn2{k}(ii), C2{k}(ii), fp2{k}(ii), W2{k}(ii)] = count(act_list, found_list2, 0);
            [tn3{k}(ii), fn3{k}(ii), C3{k}(ii), fp3{k}(ii), W3{k}(ii)] = count(act_list, found_list3, 0);
        end
        k = k + 1; % update counter for number of noise levels tested
    end

    close(h) % end waitbar
    
    for j = 1:(k-1)
        paulabars(x(j)-.025,0,tn{j},fn{j},C{j},fp{j},W{j},0);
    end
    axis([-0.025 .525 0 1]); title('Greedy with Pairs','fontsize',18)
    
    figure;
    for j = 1:(k-1)
        paulabars(x(j)-.025,0,tn2{j},fn2{j},C2{j},fp2{j},W2{j},0);
    end
    axis([-0.025 .525 0 1]); title('Simple Greedy','fontsize',18)
    
    figure;
    for j = 1:(k-1)
        paulabars(x(j)-.025,0,tn3{j},fn3{j},C3{j},fp3{j},W3{j},0);
    end
    axis([-0.025 .525 0 1]); title('Brute-Force Fitting','fontsize',18)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TIME SHIFTS %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif job == 2
    %-----------------make signal vectors--------------------------
    N = 20; % total number of sample times
    ns = 8;
    m = 4; % of spike types
    f1 = @(t, w) exp(-t.^2/(2*(w/7)^2)); % 1-peak spike
    f2 = @(t,w) 8*t.*exp(-.5*t.^2/(2*(w/7)^2)); % 2-peak spike

    for i = 1:m % pick time shifts (treat each shift as a different spike type)
        ts(1:ns,i) = 3:2:17;
    end
    for i = 1:ns 
        y(i,:) = synth(ts(i), 1, 1, 1, f1, 0, N); % spike type 1 time shifts
        if m >= 2, y(ns + i,:) = synth(ts(i), 1, 3, 1, f1, 0, N); end % spike type 2 (wider) time shifts
        if m >= 3, y(2*ns + i,:) = synth(ts(i), -.5, 1, 1, f1, 0, N); end % spike type 3 (upside down, narrower) time shifts
        if m >= 4, y(3*ns + i,:) = synth(ts(i), .8, 1, 1, f2, 0, N); end % down-up spike
        if m >= 5, y(4*ns + i,:) = synth(ts(i), -.5, 1, 1, f2, 0, N); end % up-down spike
    end
    gamma = [0.5 0.3 0.1 0.05];
    gamma2 = ones(1,ns)*gamma(1); % probability of each spike type (incl. shifts as types)
    for i = 2:m
        gamma2 = [gamma2 ones(1,ns)*gamma(i)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN EXPERIMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%
    for eta = x % signal to noise ratio
        sad = []; sad2 = []; % initialize
        for ii = 1:numii
            for i = 1:runs % runs per noise value
                % pick spike shapes present
                pick = rand(size(gamma))<gamma;
                P = 1:m;
                act_list{i} = P(pick==1); % cell array, each element = list of actual spike shapes

                % pick time shift
                ts_act{i} = rand(1,m)*N/2+5; % random uniform in [5,15]
%                 ts_act{i} = datasample(5:2:15,m); % integral time shifts

                % build signal vector 
                yn = zeros(1,N);
                for c = 1:length(act_list{i})
                    if act_list{i}(c) == 1, yn = yn + synth(ts_act{i}(c), 1, 1, 1, f1, 0, N); 
                    elseif act_list{i}(c) == 2, yn = yn + synth(ts_act{i}(c), 1, 3, 1, f1, 0, N);
                    elseif act_list{i}(c) == 3, yn = yn + synth(ts_act{i}(c), -.5, 1, 1, f1, 0, N);
                    elseif act_list{i}(c) == 4, yn = yn + synth(ts_act{i}(c), 1, 1, 1, f2, 0, N);
                    elseif act_list{i}(c) == 5, yn = yn + synth(ts_act{i}(c), -.75, 1, 1, f2, 0, N);
                    end
                end

                % add normally distributed random noise, proportional to eta
                if noise == 0, eta = 0; end
                yn = yn + eta*randn(1,N); % add noise

                % call algorithms
                if penalty == 0
                    found_list{i} = greedy_pairs(yn, y, m, gamma2, 0, repl); % cell array, each element = list of detected spike shapes
                    found_list2{i} = greedy_likely(yn, y, ns, 0, gamma2, repl);
                    found_list3{i} = brute_force_ts(yn, y, m, gamma2, 0);
                elseif penalty == 1
                    found_list{i} = greedy_pairs(yn, y, m, gamma2, eta, repl); % cell array, each element = list of detected spike shapes
                    found_list2{i} = greedy_likely(yn, y, ns, eta, gamma2, repl);
                    found_list3{i} = brute_force_ts(yn, y, m, gamma2, eta);
                end
               
                l = l + 1; % update counter for waitbar
                waitbar(l/(num_eta*runs*numii)) % update waitbar
            end
            [tn{k}(ii), fn{k}(ii), C{k}(ii), fp{k}(ii), W{k}(ii), C_ts{k}(ii), sad_temp] = count_ts2(act_list, found_list, ts, ts_act,0); % calculate error fractions
            [tn2{k}(ii), fn2{k}(ii), C2{k}(ii), fp2{k}(ii), W2{k}(ii), C_ts2{k}(ii), sad2_temp] = count_ts2(act_list, found_list2, ts, ts_act,0); % calculate error fractions
            [tn3{k}(ii), fn3{k}(ii), C3{k}(ii), fp3{k}(ii), W3{k}(ii), C_ts3{k}(ii), sad3] = count_ts2(act_list, found_list3, ts, ts_act,0); % calculate error fractions
        end

        k = k + 1; % update counter for number of noise levels tested
    end

    close(h) % end waitbar

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:(k-1)
        paulabars(x(j)-.025,0,tn{j},fn{j},C{j},fp{j},W{j},0);
    end
    axis([-0.025 .525 0 1]); title('Greedy with Pairs','fontsize',18)
    
    figure;
    for j = 1:(k-1)
        paulabars(x(j)-.025,0,tn2{j},fn2{j},C2{j},fp2{j},W2{j},0);
    end
    axis([-0.025 .525 0 1]); title('Simple Greedy','fontsize',18)
    
    figure;
    for j = 1:(k-1)
        paulabars(x(j)-.025,0,tn3{j},fn3{j},C3{j},fp3{j},W3{j},0);
    end
    axis([-0.025 .525 0 1]); title('Brute-Force Fitting','fontsize',18)
end