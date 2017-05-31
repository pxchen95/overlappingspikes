function td = detect(y, SNR, plot_on, t)
% Approximates firing times by convolving the signal vector with a
% tophat function, threshold filtering the signal, and taking the average 
% sample of each spike as the detected firing time
% 
% INPUTS: 
%     y                       = 1 x N vector, signal strength at each sample times
%     SNR                     = signal to noise ratio
%     plot_on (opt)           = 1 to show debugging plots, default is 0
%     t (req. if plot_on = 1) = 1 x L vector, actual firing times, used for plotting only
%         
% OUTPUTS: 
%     td                      = 1 x l vector, detected firing times

if nargin < 4, plot_on = 0; end % default = no plots

threshold = 6; % currently arbitrarily set threshold, goal = related to SNR
td_raw = conv(y, 3*ones(1,20), 'same'); % convolve with tophat function
td_filt = td_raw; % data filtered using threshold 
td_filt(td_filt < threshold) = 0; % signals below threshhold = 0

% estimate firing times as mean sample time of each spike
j = 1; k = 1; % initialize counters
N = length(y); % total number of samples
maxn = N; td = zeros(1,maxn); % initialize so td is nonempty
for i = 1:N
    if td_filt(i) ~= 0
        spike(j) = i; % records sample times with nonzero filtered signal
        j = j + 1;
    elseif i > 1 && td_filt(i-1) ~= 0
        td(k) = mean(spike); % detected firing times
        spike(1:end) = [ ]; % reset spike vector
        j = 1; k = k + 1;
    end
end

Nfound = k - 1; % number of detected firing events
td = td(1:Nfound); % truncate vector td
% if isempty(td), display('no spikes detected'); end % display message if no spikes detected

% plots for debugging, created if plot_switch = 1
if plot_on == 1
    L = length(t); % number of actual firing events
    subplot(3,1,1) % plot raw signal with actual firing times
    plot(1:N, y, '.-') % plot raw signal (blue)
    xlabel('sample time'); ylabel('signal strength'); title('actual signal')
    hold on
    for j = 1:L % plot green vertical line at each actual firing time
        plot([t(j),t(j)],[-1,1],'-g')
    end
    xlim([1 N])
    hold off

    subplot(3,1,2) % plot convolved signal, actual firing times, and threshold
    plot(1:N, td_raw, 'r') % plot convolved signal (red)
    xlabel('sample time'); ylabel('signal strength'); title('convolved signal')
    hold on
    for j = 1:L % plot green vertical line at each actual firing time
        plot([t(j),t(j)],[-2,20],'-g')
    end
    plot(1:N, threshold*ones(1,N)) % plot threshold (blue)
    xlim([1 N])
    hold off

    subplot(3,1,3) % plot filtered signal, detected firing times, and actual firing times
    plot(1:N, td_filt, 'r') % plot filtered signal (red)
    xlabel('sample time'); ylabel('signal strength'); title('filtered signal')
    hold on
    for j = 1:L % plot green dotted line at actual firing times
        plot(t(j)*ones(length(-2:3:30)),-2:3:30,'.g')
    end
    for j = 1:Nfound % plot blue vertical line at detected firing times
        plot([td(j),td(j)],[-2,30],'-b')
    end
    xlim([1 N])
    hold off
end

end