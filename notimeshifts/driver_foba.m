% tests foward-backward greedy algorithm

clear

repl = 0; % 1 = greedy with replacement, 0 = greedy without replacement
noise = 0; % 0 = no noise, 1 = with noise
use_bars = 0; % 0 = make scatter plot for no noise, 1 = make bar graph for no noise

%-----------------make signal vectors--------------------------
N = 20; % total number of sample times
f1 = @(t, w) exp(-t.^2/(2*(w/7)^2)); % 1-peak spike
y(1,:) = synth(10, 1, 1, 1, f1, 0, N); % spike type 1
y(2,:) = synth(10, -0.8, 1, 1, f1, 0, N); % spike type 2 (wider spike)
y(3,:) = (y(1,:) + y(2,:))*1.1; % upside down, short spike
[m,n] = size(y); % # of spike types vs # of samples
gamma = ones(1,m)*0.5; % probability of each spike type
%--------------------------------------------------------------

x = 0:0.1:10; % range of noise levels tested 
num_eta = length(x); % number of noise levels tested
k = 1; % counter for number of noise levels tested
h = waitbar(0,'Please wait...'); % make wait bar

for eta = x % signal to noise ratio
    for i = 1:1000 % 1000 runs per noise value     
        % pick spike shapes present
        pick = rand(size(gamma))<gamma;
        P = 1:m;
        act_list{i} = P(pick==1); % cell array, each element = list of actual spike shapes
        
        % add normally distributed random noise, proportional to eta
        if noise == 0, eta = 0; end
        
        if sum(pick) > 1 % create signal vector
            yn = sum(y(act_list{i},:)) + eta*randn(1,N); 
        elseif sum(pick) == 0
            yn = eta*randn(1,N);
        else
            yn = y(act_list{i},:) + eta*randn(1,N);
        end
        
        brute_found{i} = brute_force(yn, y);
        foba_found{i} = FoBa(yn, y, repl);
        forward_found{i} = greedy_sort(yn, y, repl);
    end
    
    [tnb(k), fnb(k), Cb(k), fpb(k), Wb(k)] = count(act_list, brute_found);
    [tnf(k), fnf(k), Cf(k), fpf(k), Wf(k)] = count(act_list, foba_found);
    [tn(k), fn(k), C(k), fp(k), W(k)] = count(act_list, forward_found);
    
    k = k + 1;
    waitbar(k/num_eta) % update waitbar
end

close(h) % end waitbar 

% plot error fractions
plot(x,tnb,'--b', x,fnb,'--c', x,Cb,'--g', x,fpb,'--m', x,Wb,'--r')
hold on
axis tight; ylim([0 1])
if noise == 0, xlabel('# of runs'); 
else xlabel('noise (eta = 1/SNR)'); end
ylabel('error fraction')
legend('true negative', 'false negative', 'correct', ...
    'false positive', 'wrong', 'location', 'NorthEast')
plot(x,tnf,'.b', x,fnf,'.c', x,Cf,'.g', x,fpf,'.m', x,Wf,'.r')
plot(x,tn,'b', x,fn,'c', x,C,'g', x,fp,'m', x,W,'r')
title('-- brute force, \cdot FoBa, - simple greedy', 'fontsize', 18)
hold off