% runs and tests synth.m

clear 

N = 200; % total number of sample times
L = 10; % desired number of events
refrac = 20; % refractory period

%rng(1) % fixes seed
t = gen_firing(N,L,refrac); % call gen_firing.m to make firing times
display('firing times'); display(t)
L = length(t); % actual number of firing events
display('# of firing events'); display(L)

a = ones(1,L); % amplitudes = 1
w = ones(1,L); % spike widths = "1"
s = ones(1,L);  % spike types = 1

f1 = @(t, w) exp(-t.^2/(2*(w/7)^2)); % 1-peak spike
f2 = @(t, w) exp(-t.^2/(2*(w/7)^2)).*t/(exp(-1/2)*(w/7)); % 2-peak spike

% function tests
% x = -1:0.1:1;
% plot(x, f1(x,2));
% plot(x, f2(x,2));

y = synth(t, a, w, s, f1, f2, N); % call synth.m to create signal vector
snr = 3; % signal-to-noise ratio
eta = 0.3;%1/snr; % noise level (SNR = 1/eta)
y = y + eta*randn(1,N); % add normally distributed random noise, proportional to eta

% plot sample time vs signal strength
figure
plot(1:N, y, '.-','linewidth',1.5,'markersize',15) 
xlabel('Sample Time', 'fontsize', 18); ylabel('Signal Amplitude', 'fontsize', 18)
set(gca,'fontsize',18)

% check that spikes are placed at the correct firing times
hold on
for j = 1:L % plot green vertical line at each t(j)
    plot([t(j),t(j)],[-1,1],'-g','linewidth',1.5);
end
xlim([1 N])
legend('Synthesized Signal','Firing Time')
hold off