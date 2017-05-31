% sort synthetic clips using templates spikes from harrisdata (with and 
% without time shifts)

clear

job = 2; % 1 = without time shifts, 2 = with time shifts
noise = 1; repl = 0;
rng(0); % fix seed 

%%%%%%%%%%%%%% LOAD DATA %%%%%%%%%%%%%%
file = matfile('harrisdata.mat');
K = file.K; % # of spike types
Y = file.Y; % filtered signal
clips = file.clips; % clips to be sorted
clipsize = file.clipsize; % sample times per clip
firingtimes = file.firingtimes;
groundtruthlabels = file.groundtruthlabels; % which clips contain spike 2
nclips = file.nclips; % # of clips to be sorted
ntrue = file.ntrue;
samples_per_sec = file.samples_per_sec;
sortedlabels = file.sortedlabels; % spikes detected by different algorithm
templates = file.templates; % known spike types
truetimes = file.truetimes;

% --------------- applies for all jobs --------------- %

x = 0:0.1:0.5; % range of noise levels tested
runs = 200;    
numii = 5; %50;     % # of runs per noise level = runs * numii
num_eta = length(x); % number of noise levels tested
k = 1; % counter for number of noise levels tested
l = 1; % counter for waitbar
j = 1; % counter for total number of runs
h = waitbar(0,'Please wait...'); % make wait bar
SAD = []; SAD2 = []; SAD3 = []; % initialize, standard absolute deviation

% ---------------------- begin cases ----------------- %

y = templates;
ymax = max(abs(y(:))); y = y/ymax;
m = K; N = 30;
gamma = ones(1,m)*0.5;
if job == 1
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
                found_list{i} = greedy_pairs(yn, y, m, gamma, eta, repl); % cell array, each element = list of detected spike shapes
                found_list2{i} = greedy_likely(yn, y, 1, eta, gamma, repl);
                found_list3{i} = brute_force_ts(yn, y, m, gamma, eta);

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
    N = clipsize; % total number of sample times
    temp = 9:2:21; % time shifts tested
    ns = length(temp); % number of times shifts tested
    ts = []; % initialize
    for i = 1:K
        ts = [ts temp]; % time shifts to test per spike type
    end

    m = length(temp); % total # of spike types tested (incl. shifts as types)
    ftemp = templates/max(abs(templates(:)));
    fmax = max(abs(templates(:)));
    f = zeros(m,clipsize); 
    for i = 1:m
        for j = 1:K
            if ts(i) < 15 % time shift to the left
                last = clipsize - (15 - ts(i)) + 1;
                f((j-1)*m + i,1:last) = ftemp(j,(15-ts(i)):clipsize);
            elseif ts(i) == 15 % same as template ("centered")
                f((j-1)*m + i,:) = ftemp(j,:);
            else % time shift to the right
                first = 1 + ts(i) - 15;
                f((j-1)*m + i,first:clipsize) = ftemp(j,1:(clipsize - first)+1);
            end
        end
    end
    m = K; % of spike types
    gamma = 0.5*ones(1,3); %[0.1 0.25 0.5];
    gamma2 = ones(1,ns)*gamma(1); % probability of each spike type (incl. shifts as types)
    for i = 2:m
        gamma2 = [gamma2 ones(1,ns)*gamma(i)];
    end
    
    ts = zeros(ns,m);
    for i = 1:m % pick time shifts (treat each shift as a different spike type)
        ts(1:ns,i) = temp;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN EXPERIMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%
    G = 1; 
    for eta = x % signal to noise ratio
        for ii = 1:numii
            for i = 1:runs % runs per noise value
                % pick spike shapes present
                pick = rand(size(gamma))<gamma;
                P = 1:m;
                act_list{i} = P(pick==1); % cell array, each element = list of actual spike shapes

                % pick time shift
                ts_act{i} = datasample(11:1:19,m);

                % build signal vector 
                yn = zeros(1,N);
                for c = 1:length(act_list{i})
                    if ts_act{i}(act_list{i}(c)) < 15 % time shift to the left
                        last = clipsize - (15 - ts_act{i}(act_list{i}(c))) + 1;
                        yn(1,1:last) = yn(1,1:last) + ftemp(act_list{i}(c),(15-ts_act{i}(act_list{i}(c))):clipsize);
                    elseif ts_act{i}(act_list{i}(c)) == 15 % same as template ("centered")
                        yn = yn + ftemp(act_list{i}(c),:);
                    else % time shift to the right
                        first = 1 + ts_act{i}(act_list{i}(c)) - 15;
                        yn(1,first:clipsize) = yn(1,first:clipsize) + ftemp(act_list{i}(c),1:(clipsize - first)+1);
                    end
                end

                % add normally distributed random noise, proportional to eta
                if noise == 0, eta = 0; end
                yn = yn + eta*randn(1,N); % add noise

                % forward greedy algorithm
                found_list{i} = greedy_pairs(yn, f, m, gamma2, eta, repl); % cell array, each element = list of detected spike shapes
                found_list2{i} = greedy_likely(yn, f, ns, eta, gamma2, repl); % cell array, each element = list of detected spike shapes
                found_list3{i} = brute_force_ts(yn, f, m, gamma2, eta);               

                l = l + 1; % update counter for waitbar
                waitbar(l/(num_eta*runs*numii)) % update waitbar
            end
            [tn{k}(ii), fn{k}(ii), C{k}(ii), fp{k}(ii), W{k}(ii), C_ts{k}(ii), sad_temp] = count_ts2(act_list, found_list, ts, ts_act,0); % calculate error fractions
            [tn2{k}(ii), fn2{k}(ii), C2{k}(ii), fp2{k}(ii), W2{k}(ii), C_ts2{k}(ii), sad2_temp] = count_ts2(act_list, found_list2, ts, ts_act,0); % calculate error fractions
            [tn3{k}(ii), fn3{k}(ii), C3{k}(ii), fp3{k}(ii), W3{k}(ii), C_ts3{k}(ii), sad3] = count_ts2(act_list, found_list3, ts, ts_act,0); % calculate error fractions

            % for looking at clips for which the algorithms have sorting errors
%             if fp3{k}(ii) > fp2{k}(ii)
%                 eta
%                 ywn{G} = yn; 
%                 act_list{1}
%                 ts_act{1}
%                 ceil(found_list{1}/ns)
% %                 found_list{1}
%                 ceil(found_list2{1}/ns)
%                 found_list2{1}
%                 ceil(found_list3{1}/ns)
%                 found_list3{1}
%                 G = G + 1;
%             end
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