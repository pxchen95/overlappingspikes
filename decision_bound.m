% plots decision boundaries for TWO spike case for the greedy algorithm
% (without replacement) and the brute-force algorithm

clear

noise = 0; % 0 = no noise, 1 = with noise

%-----------------make signal vectors--------------------------
N = 20; % total number of sample times
f1 = @(t, w) exp(-t.^2/(2*(w/7)^2)); % 1-peak spike
% f2 = @(t, w) t.*exp(-t.^2/(2*(w/7)^2)); % 2-peak spike
y(1,:) = synth(10, 1, 1, 1, f1, 0, N); % spike type 1
% y(2,:) = synth(10, -.5, 0.5, 1, f1, 0, N); % upside down, short spike
y(2,:) = synth(10, 1, 3, 1, f1, 0, N); % spike type 2 (wider spike)
y(3,:) = y(1,:) + y(2,:); % spike type 1 + spike type 2
%--------------------------------------------------------------

range = 0:0.2:5; % range of noise levels tested 
num_eta = length(range); % number of noise levels tested
k = 1; % counter for number of noise levels tested
num_W = 0; % counts number of wrong identifications
total = 0; % counts total number of detected spikes
h = waitbar(0,'Please wait...'); % make wait bar

% orthonormalize vectors via Gram-Schmidt process
u1 = y(1,:); e1 = u1/norm(u1); 
u2 = y(2,:) - dot(y(2,:),u1)/dot(u1,u1) * u1; e2 = u2/norm(u2); 

% project no noise spike types into 2D plane
x = [0, dot(y(1,:),e1), dot(y(2,:),e1), dot(y(3,:),e1)];
yy = [0, dot(y(1,:),e2), dot(y(2,:),e2), dot(y(3,:),e2)];

% begin plotting 
% subplot(2,3,1) % plots actual spikes
% plot(x,yy,'k*') % plot ("actual spike, no noise")
% title('Actual Spikes')
% hold on

subplot(1,2,1) % plots detection by brute-force
plot(x,yy,'k.') % plot ("actual spike, no noise")
title('Decision Boundaries for Brute-Force Fitting')
hold on

subplot(1,2,2) % plots detection by foward greedy
plot(x,yy,'k.') % plot ("actual spike, no noise")
title('Decision Boundaries for Simple Forward Greedy')
hold on

% subplot(2,3,5) % plots detection by backward greedy
% plot(x,yy,'k*') % plot ("actual spike, no noise")
% title('Decision Boundaries for Backward Greedy Algorithm')
% hold on

% subplot(2,3,6) % plots detection by FoBa greedy
% plot(x,yy,'k*') % plot ("actual spike, no noise")
% title('Decision Boundaries for FoBa Greedy Algorithm')
% hold on

for eta = range % signal to noise ratio
    for i = 1:100 % runs per noise value
        
        % actual spikes
        y_0 = eta*randn(1,N); % no spike
        y_1 = y(1,:) + eta*randn(1,N); % spike type 1      
        y_2 = y(2,:) + eta*randn(1,N); % spike type 2
        y_3 = y(3,:) + eta*randn(1,N); % spike types 1 & 2
        
        % detect spikes via foward greedy algorithm (no replacement)
        y0_det = greedy_sort(y_0, y(1:2,:), 0);
        y1_det = greedy_sort(y_1, y(1:2,:), 0);
        y2_det = greedy_sort(y_2, y(1:2,:), 0);
        y3_det = greedy_sort(y_3, y(1:2,:), 0);
        
        % detect spikes via backward greedy algorithm
%         y0_det = greedy_back(y_0, y(1:2,:));
%         y1_det = greedy_back(y_1, y(1:2,:));
%         y2_det = greedy_back(y_2, y(1:2,:));
%         y3_det = greedy_back(y_3, y(1:2,:));
        
        % detect spikes via brute-force algorithm
        y0_detb = brute_force(y_0, y(1:2,:));
        y1_detb = brute_force(y_1, y(1:2,:));
        y2_detb = brute_force(y_2, y(1:2,:));
        y3_detb = brute_force(y_3, y(1:2,:));
%         
%         % detect spikes via FoBa greedy algorithm
%         y0_detf = FoBa(y_0, y(1:2,:), 0);
%         y1_detf = FoBa(y_1, y(1:2,:), 0);
%         y2_detf = FoBa(y_2, y(1:2,:), 0);
%         y3_detf = FoBa(y_3, y(1:2,:), 0);        
        
        % plotting actual spikes
%         subplot(2,3,1)
%         plot(dot(y_0,e1),dot(y_0,e2),'.b', dot(y_1,e1),dot(y_1,e2),'.g', ...
%             dot(y_2,e1),dot(y_2,e2),'.r', dot(y_3,e1),dot(y_3,e2),'.m')
        
        % plotting brute-force
        subplot(1,2,1) 
        if isempty(y0_detb), plot(dot(y_0,e1),dot(y_0,e2),'.b');
        elseif y0_detb == 1, plot(dot(y_0,e1),dot(y_0,e2),'.g');
        elseif y0_detb == 2, plot(dot(y_0,e1),dot(y_0,e2),'.r');
        else plot(dot(y_0,e1), dot(y_0,e2),'.m');
        end
        
        if isempty(y1_detb), plot(dot(y_1,e1),dot(y_1,e2),'.b');
        elseif y1_detb == 1, plot(dot(y_1,e1),dot(y_1,e2),'.g');
        elseif y1_detb == 2, plot(dot(y_1,e1),dot(y_1,e2),'.r');
        else plot(dot(y_1,e1), dot(y_1,e2),'.m');
        end
        
        if isempty(y2_detb), plot(dot(y_2,e1),dot(y_2,e2),'.b');
        elseif y2_detb == 1, plot(dot(y_2,e1),dot(y_2,e2),'.g'); 
        elseif y2_detb == 2, plot(dot(y_2,e1),dot(y_2,e2),'.r'); 
        else plot(dot(y_2,e1), dot(y_2,e2),'.m'); 
        end
        
        if isempty(y3_detb), plot(dot(y_3,e1),dot(y_3,e2),'.b');
        elseif y3_detb == 1, plot(dot(y_3,e1),dot(y_3,e2),'.g');
        elseif y3_detb == 2, plot(dot(y_3,e1),dot(y_3,e2),'.r'); 
        else plot(dot(y_3,e1), dot(y_3,e2),'.m');
        end
        
        % plotting foward greedy
        subplot(1,2,2)
        if isempty(y0_det), plot(dot(y_0,e1),dot(y_0,e2),'.b');
        elseif y0_det == 1, plot(dot(y_0,e1),dot(y_0,e2),'.g');
        elseif y0_det == 2, plot(dot(y_0,e1),dot(y_0,e2),'.r'); 
        else plot(dot(y_0,e1), dot(y_0,e2),'.m');
        end
        
        if isempty(y1_det), plot(dot(y_1,e1),dot(y_1,e2),'.b');
        elseif y1_det == 1, plot(dot(y_1,e1),dot(y_1,e2),'.g');
        elseif y1_det == 2, plot(dot(y_1,e1),dot(y_1,e2),'.r');
        else plot(dot(y_1,e1), dot(y_1,e2),'.m');
        end
        
        if isempty(y2_det), plot(dot(y_2,e1),dot(y_2,e2),'.b');
        elseif y2_det == 1, plot(dot(y_2,e1),dot(y_2,e2),'.g'); 
        elseif y2_det == 2, plot(dot(y_2,e1),dot(y_2,e2),'.r'); 
        else plot(dot(y_2,e1), dot(y_2,e2),'.m');
        end
        
        if isempty(y3_det), plot(dot(y_3,e1),dot(y_3,e2),'.b');
        elseif y3_det == 1, plot(dot(y_3,e1),dot(y_3,e2),'.g');
        elseif y3_det == 2, plot(dot(y_3,e1),dot(y_3,e2),'.r'); 
        else plot(dot(y_3,e1), dot(y_3,e2),'.m'); 
        end
        
        % plotting backward greedy
%         subplot(2,3,5)
%         if isempty(y0_det2), plot(dot(y_0,e1),dot(y_0,e2),'.b');
%         elseif y0_det2 == 1, plot(dot(y_0,e1),dot(y_0,e2),'.g');
%         elseif y0_det2 == 2, plot(dot(y_0,e1),dot(y_0,e2),'.r'); 
%         else plot(dot(y_0,e1),dot(y_0,e2),'.m'); 
%         end
%         
%         if isempty(y1_det2), plot(dot(y_1,e1),dot(y_1,e2),'.b');
%         elseif y1_det2 == 1, plot(dot(y_1,e1),dot(y_1,e2),'.g'); 
%         elseif y1_det2 == 2, plot(dot(y_1,e1),dot(y_1,e2),'.r');
%         else plot(dot(y_1,e1),dot(y_1,e2),'.m'); 
%         end
%         
%         if isempty(y2_det2), plot(dot(y_2,e1),dot(y_2,e2),'.b');
%         elseif y2_det2 == 1, plot(dot(y_2,e1),dot(y_2,e2),'.g');
%         elseif y2_det2 == 2, plot(dot(y_2,e1),dot(y_2,e2),'.r'); 
%         else plot(dot(y_2,e1),dot(y_2,e2),'.m');
%         end
%         
%         if isempty(y3_det2), plot(dot(y_3,e1),dot(y_3,e2),'.b');
%         elseif y3_det2 == 1, plot(dot(y_3,e1),dot(y_3,e2),'.g'); 
%         elseif y3_det2 == 2, plot(dot(y_3,e1),dot(y_3,e2),'.r');
%         else plot(dot(y_3,e1),dot(y_3,e2),'.m');
%         end
        
        % plotting FoBa
%         subplot(2,3,6) 
%         if isempty(y0_detf), plot(dot(y_0,e1),dot(y_0,e2),'.b');
%         elseif y0_detf == 1, plot(dot(y_0,e1),dot(y_0,e2),'.g');
%         elseif y0_detf == 2, plot(dot(y_0,e1),dot(y_0,e2),'.r');
%         else plot(dot(y_0,e1), dot(y_0,e2),'.m');
%         end
%         
%         if isempty(y1_detf), plot(dot(y_1,e1),dot(y_1,e2),'.b');
%         elseif y1_detf == 1, plot(dot(y_1,e1),dot(y_1,e2),'.g');
%         elseif y1_detf == 2, plot(dot(y_1,e1),dot(y_1,e2),'.r');
%         else plot(dot(y_1,e1), dot(y_1,e2),'.m');
%         end
%         
%         if isempty(y2_detf), plot(dot(y_2,e1),dot(y_2,e2),'.b');
%         elseif y2_detf == 1, plot(dot(y_2,e1),dot(y_2,e2),'.g'); 
%         elseif y2_detf == 2, plot(dot(y_2,e1),dot(y_2,e2),'.r'); 
%         else plot(dot(y_2,e1), dot(y_2,e2),'.m'); 
%         end
%         
%         if isempty(y3_detf), plot(dot(y_3,e1),dot(y_3,e2),'.b');
%         elseif y3_detf == 1, plot(dot(y_3,e1),dot(y_3,e2),'.g');
%         elseif y3_detf == 2, plot(dot(y_3,e1),dot(y_3,e2),'.r'); 
%         else plot(dot(y_3,e1), dot(y_3,e2),'.m');
%         end
    end
    
    k = k + 1; % update counter
    waitbar(k/num_eta) % update waitbar
end

close(h) % close waitbar

% subplot(2,3,1)
% plot(x,yy,'k*') % plot ("actual spike, no noise")
% legend('No Noise','No Spike','Spike Type 1','Spike Type 2',...
%     'Spike Type 1 + 2','Location','northwest')
% hold off

subplot(1,2,1)
plot(x,yy,'k.','MarkerSize', 24) % plot ("actual spike, no noise")
legend('Actual Types (No Noise)','Detect No Spike','Detect Type 1','Detect Type 3',...
     'Detect Types 1 + 3','Location','northwest')
xlim([-5,5]); ylim([-5,5])
% set(gca,'fontsize',36)
hold off

subplot(1,2,2)
plot(x,yy,'k.','MarkerSize', 24) % plot ("actual spike, no noise")
% legend('No Noise','No Spike','Spike Type 1','Spike Type 2',...
%     'Spike Type 1 + 2','Location','northwest')
xlim([-5,5]); ylim([-5,5])
% set(gca,'fontsize',36)
hold off

% subplot(2,3,5)
% plot(x,yy,'k*') % plot ("actual spike, no noise")
% % legend('No Noise','No Spike','Spike Type 1','Spike Type 2',...
% %     'Spike Type 1 + 2','Location','northwest')
% hold off
% 
% subplot(2,3,6)
% plot(x,yy,'k*') % plot ("actual spike, no noise")
% % legend('No Noise','No Spike','Spike Type 1','Spike Type 2',...
% %     'Spike Type 1 + 2','Location','northwest')
% hold off

% code for centering plots with odd number of subplots
% z(1) = subplot(2,2,1); z(2) = subplot(2,2,3); z(3) = subplot(2,2,4); 
% pos = get(z,'Position');
% new = mean(cellfun(@(v)v(1),pos(2:3)));
% set(z(1),'Position',[new,pos{1}(2:end)]);

% plot spike types
% plot(y(1,:),'b'); hold on; plot(y(2,:),'g'); hold off