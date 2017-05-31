% run algorithms on harris data (real clips)
clear 

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

%%%%%%%%%%%%%% MAKE TIME SHIFTS TO TEST %%%%%%%%%%%%%%
temp = 9:2:21; % number of times shifts tested
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

m = length(ts);

%%%%%%%%%%%%%% RUN ALGORITHMS 1 %%%%%%%%%%%%%%
eta = 0.2; repl = 0; gamma = 0.5*ones(1,m); %gamma((length(temp)+1):(2*length(temp))) = 0.2790;
h = waitbar(0,'Step 1/3: Spike Sorting...'); % make wait bar
y = transpose(clips); % each row = 1 clip
y = y/fmax;
for i = 1:nclips
        found_list1{i} = brute_force_ts(y(i,:), f, K, gamma, eta);
        found_list3{i} = greedy_likely(y(i,:), f, m, eta, gamma, repl);
        found_list2{i} = greedy_pairs(y(i,:), f, K, gamma, eta, repl);
        waitbar(i/nclips) % update waitbar
end
close(h) % end waitbar

%%%%%%%%%%%%%% ERROR CALCULATION 1 %%%%%%%%%%%%%%
h = waitbar(0,'Step 2/3: Error Analysis 1...'); % make wait bar
ts_tested = vertcat(temp, temp); ts_tested = vertcat(ts_tested, temp);
ts_tested = transpose(ts_tested);
ts_act = 15*ones(1,nclips);
tn1 = zeros(1,5); fn1 = zeros(1,5); C1 = zeros(1,5); fp1 = zeros(1,5); W1 = zeros(1,5);
tn2 = zeros(1,5); fn2 = zeros(1,5); C2 = zeros(1,5); fp2 = zeros(1,5); W2 = zeros(1,5);
tn3 = zeros(1,5); fn3 = zeros(1,5); C3 = zeros(1,5); fp3 = zeros(1,5); W3 = zeros(1,5);
for k = 1:5
    for i = 1:529
        [tn, fn, C, fp, W] = count_ts2(sortedlabels((k-1)*529 + i), found_list1{(k-1)*529 + i}, ts_tested, ts_act, 0); 
        tn1(k) = tn1(k) + tn; fn1(k) = fn1(k) + fn; C1(k) = C1(k) + C; 
        fp1(k) = fp1(k) + fp; W1(k) = W1(k) + W;
        
        [tn, fn, C, fp, W] = count_ts2(sortedlabels((k-1)*529 + i), found_list2{(k-1)*529 + i}, ts_tested, ts_act, 0); 
        tn2(k) = tn2(k) + tn; fn2(k) = fn2(k) + fn; C2(k) = C2(k) + C; 
        fp2(k) = fp2(k) + fp; W2(k) = W2(k) + W;
        
        [tn, fn, C, fp, W] = count_ts2(sortedlabels((k-1)*529 + i), found_list3{(k-1)*529 + i}, ts_tested, ts_act, 0); 
        tn3(k) = tn3(k) + tn; fn3(k) = fn3(k) + fn; C3(k) = C3(k) + C; 
        fp3(k) = fp3(k) + fp; W3(k) = W3(k) + W;
        
        waitbar((((k-1)*529)+i)/nclips) % update waitbar
    end
end
close(h) % end waitbar

%%%%%%%%%%%%%% PLOTTING ERROR BARS 1 %%%%%%%%%%%%%%
figure
x = 1:3;
tnk{1} = tn1; tnk{2} = tn2; tnk{3} = tn3;
fnk{1} = fn1; fnk{2} = fn2; fnk{3} = fn3;
Ck{1} = C1; Ck{2} = C2; Ck{3} = C3;
fpk{1} = fp1; fpk{2} = fp2; fpk{3} = fp3;
Wk{1} = W1; Wk{2} = W2; Wk{3} = W3;
for j = 1:3
    paulabars2(x(j)-.25,0,tnk{j},fnk{j},Ck{j},fpk{j},Wk{j},0);
end
axis([0.5 3.5 0 1])

%%%%%%%%%%%%%% REDEFINE ACT_LIST/FOUND_LISTS %%%%%%%%%%%%%%
ns = length(temp); 
for i = 1:nclips
    if isnan(groundtruthlabels(i)), act_list{i} = [];
    else act_list{i} = groundtruthlabels(i);
    end
    
    if ~isempty(found_list1{i})
        if sum(ceil(found_list1{i}/ns) == 2) == 0, found_list1{i} = []; % only count if 2 was detected
        elseif length(found_list1{i}) >= 1, found_list1{i} = 2; % neglect any extra spikes found
        end
    end
    
    if ~isempty(found_list2{i})
        if sum(ceil(found_list2{i}/ns) == 2) == 0, found_list2{i} = [];
        elseif length(found_list2{i}) >= 1, found_list2{i} = 2;
        end
    end
    
    if ~isempty(found_list3{i})
        if sum(ceil(found_list3{i}/ns) == 2) == 0, found_list3{i} = [];
        elseif length(found_list3{i}) >= 1, found_list3{i} = 2;
        end
    end
end

%%%%%%%%%%%%%% ERROR CALCULATION 2 %%%%%%%%%%%%%%
h = waitbar(0,'Step 3/3: Error Analysis 2...'); % make wait bar
ts_tested = vertcat(temp, temp); ts_tested = vertcat(ts_tested, temp);
ts_tested = transpose(ts_tested);
ts_act = 15*ones(1,nclips);
tn1 = zeros(1,5); fn1 = zeros(1,5); C1 = zeros(1,5); fp1 = zeros(1,5); W1 = zeros(1,5);
tn2 = zeros(1,5); fn2 = zeros(1,5); C2 = zeros(1,5); fp2 = zeros(1,5); W2 = zeros(1,5);
tn3 = zeros(1,5); fn3 = zeros(1,5); C3 = zeros(1,5); fp3 = zeros(1,5); W3 = zeros(1,5);
for k = 1:5
    for i = 1:529
        [tn, fn, C, fp, W] = count(act_list{(k-1)*529 + i}, found_list1{(k-1)*529 + i}, 0);
        tn1(k) = tn1(k) + tn; fn1(k) = fn1(k) + fn; C1(k) = C1(k) + C; 
        fp1(k) = fp1(k) + fp; W1(k) = W1(k) + W;
        
        [tn, fn, C, fp, W] = count(act_list{(k-1)*529 + i}, found_list2{(k-1)*529 + i}, 0);
        tn2(k) = tn2(k) + tn; fn2(k) = fn2(k) + fn; C2(k) = C2(k) + C; 
        fp2(k) = fp2(k) + fp; W2(k) = W2(k) + W;
        
        [tn, fn, C, fp, W] = count(act_list{(k-1)*529 + i}, found_list3{(k-1)*529 + i}, 0);
        tn3(k) = tn3(k) + tn; fn3(k) = fn3(k) + fn; C3(k) = C3(k) + C; 
        fp3(k) = fp3(k) + fp; W3(k) = W3(k) + W;
        
        waitbar((((k-1)*529)+i)/nclips) % update waitbar
    end
end
close(h) % end waitbar

%%%%%%%%%%%%%% PLOTTING ERROR BARS 2 %%%%%%%%%%%%%%
figure
x = 1:3;
tnk{1} = tn1; tnk{2} = tn2; tnk{3} = tn3;
fnk{1} = fn1; fnk{2} = fn2; fnk{3} = fn3;
Ck{1} = C1; Ck{2} = C2; Ck{3} = C3;
fpk{1} = fp1; fpk{2} = fp2; fpk{3} = fp3;
Wk{1} = W1; Wk{2} = W2; Wk{3} = W3;
for j = 1:3
    paulabars2(x(j)-.25,0,tnk{j},fnk{j},Ck{j},fpk{j},Wk{j},0);
end
axis([0.5 3.5 0 1])