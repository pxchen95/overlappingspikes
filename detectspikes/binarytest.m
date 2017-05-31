% Runs binary test - detect when either one or no spike is present
% plots:
%        (1) error fractions vs SNRs
%        (2) diagram of optimal threshold and pdfs for the detection parameter

clear 

%-----------------make signal vectors--------------------------
N = 20; % total number of sample times
f1 = @(t, w) exp(-t.^2/(2*(w/7)^2)); % 1-peak spike
y1 = synth(10, 1, 1, 1, f1, 0, N); % call synth.m to create 1 spike signal vector
y0 = zeros(1,N); % 0 spike signal vector
%--------------------------------------------------------------

x = 0.1:0.1:5; % range of SNRs tested
num_SNR = length(x); % number of SNRs tested
k = 1; % counter, counts number of SNRs tested
H = 1/2*dot(y1,y1); % threshold, signal above this = spike detected
h = waitbar(0,'Please wait...'); % make wait bar

for SNR = x % signal to noise ratio
    count = zeros(2,2); % 2 x 2 truth x detection matrix, [TN, FP; FN, TP]

    for i = 1:2000 % 1000 runs per SNR
        is_spike = round(rand); % randomly pick: 0 = no spike, 1 = 1 spike
        if is_spike == 0, y = y0; else y = y1; end 
        
        eta = 1/SNR; % noise level (SNR = 1/eta)
        y2 = y + eta*randn(1,N); % add normally distributed random noise, proportional to eta
        
        z = dot(y2,y1); % dot signal with spike shape
        det_spike = z>H; % 0 = no spike detected, 1 = spike detected
        
        % update count matrix 
        if is_spike < det_spike % is_spike = 0, det_spike = 1 --- detect spike that is not there (FP)
            count(1,2) = count(1,2) + 1;
        elseif is_spike > det_spike % is_spike = 1, det_spike = 0 --- miss the spike (FN)
            count(2,1) = count(2,1) + 1; 
        elseif is_spike == 0 % is_spike = det_spike = 0, no spike, none detected (TN)
            count(1,1) = count(1,1) + 1;
        else % is_spike = det_spike = 1, spike is correctly detected (TP)
            count(2,2) = count(2,2) + 1; 
        end
    end
    
    % calculate error fractions via row normalization 
    if sum(count(1,:)) ~= 0, count(1,:) = count(1,:)/sum(count(1,:)); end % divide by # of non-spikes
    if sum(count(2,:)) ~= 0, count(2,:) = count(2,:)/sum(count(2,:)); end % divide by # of true spikes
    
    tn(k) = count(1,1); fp(k) = count(1,2); % record error fractions for each SNR
    fn(k) = count(2,1); tp(k) = count(2,2); 
    
    k = k + 1; % update counter
    waitbar(k/num_SNR) % update waitbar
end

close(h) % end waitbar

% plot error fractions vs SNR
figure
plot(x,tn,'.b', x,fp,'.c', x,fn,'.m', x,tp,'.g','markersize',20)
xlabel('SNR','fontsize',18); ylabel('Error Fraction','fontsize',18)
legend('True Negative', 'False Positive (Type I)', 'Missed (Type II)', 'Correct', 'location', 'east')
hold on

% plot theoretical curves on same graph
sigma = 1./x*norm(y1); % standard deviation
predict_tn = normcdf(H,0,sigma);
predict_fn = normcdf(H,norm(y1)^2,sigma);
predict_fp = 1-predict_tn;
predict_tp = 1-predict_fn;
plot(x,predict_tn,'b', x,predict_fp,'c', x,predict_fn,'m', x,predict_tp,'g','linewidth',3)
set(gca,'fontsize',18)
hold off

% print count matrix of truth x detection (currently prints final fractions)
M = [tn(end), fp(end); fn(end), tp(end)];
printmat(M, 'Count Matrix', 'Truth=0 Truth=1', 'Detected=0 Detected=1')

% plot p(z)
figure
x_z = -1.5:0.1:4;
p1 = normpdf(x_z, 0, sigma(40)); p2 = normpdf(x_z, dot(y1,y1), sigma(40));
plot(x_z,p1, x_z,p2,'linewidth',2); ylabel('pdf','fontsize',18)
text(-1.5,0.85,'p(\Psi) (no spike) \rightarrow','fontsize',18); text(2.75,0.85,'\leftarrow p(\Psi) (1 spike)','fontsize',18)
hold on
H_z = 1/2*dot(y1,y1);
plot([H_z,H_z], [0,1], '--r','linewidth',2)
text(H_z-0.12, 0.5, '\Omega_{opt}', 'rotation',90,'fontsize',18)
text(0,0.05,'\Psi < \Omega_{opt}', 'HorizontalAlignment','center','fontsize',18)
text(H_z*2,0.05,'\Psi > \Omega_{opt}', 'HorizontalAlignment','center','fontsize',18)
set(gca,'fontsize',18)
hold off
axis tight