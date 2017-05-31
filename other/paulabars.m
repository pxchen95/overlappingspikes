function h = paulabars(xo,yo,tn,fn,cor,fp,wr,tswitch)
% PAULABARS  overlays two bar-plots showing error fractions in spike sorting
%
% h = paulabars(xo,yo,tn,fn,cor,fp,wr)
% plots a bar-like object showing the fractions of errors of each type, given
% the raw counts:
%  tn = # true negative
%  fn = # false negative (missed)
%  cor = # correct (true positive, correct identity)
%  fp = # false positive
%  wr = # wrong (true positive, wrong identity)
% (xo,yo) is the offset in the plot (lets it be added at different locations).
% tswitch: 1 = with text, 0 = without text
% The size of the bar is 0.5 across by 1 tall.
%
%  If the counts are each lists, it assumes they are indep batches of runs
%  and uses it to estimate and plot 1-sigma errorbars (assuming normal).
%
% Called without arguments, it does a self-demo.
%
% Barnett 4/22/17

if nargin==0, demo_paulabars; return; end

N = numel(tn);   % number of samples (run batches)
if numel(fn)~=N || numel(cor)~=N || numel(fp)~=N || numel(wr)~=N
  error('input counts must have same number of elements!');
end

n_det = cor+fp+wr; % # detected spikes
n_act = cor+fn+wr; % # actual spikes
n_emp = tn+fp;     % # clips actually empty

% define our fractions (for each batch if N>1)...
sens = cor./n_act;  % note difference from f_cor
f_miss = fn./n_act;
f_wr = wr./n_act;
f_fp = fp./n_det;
f_wrd = wr./n_det;
f_cor = cor./n_det;

% means for plotting
m_sens = mean(sens);
m_miss = mean(f_miss);
m_wr = mean(f_wr);
m_fp = mean(f_fp);
m_wrd = mean(f_wrd);
m_cor = mean(f_cor);

% first bar...
w = .05;   % bar width (height will be 1) %%%%%%%%%%%%%%%%%%%%%% **************
% sens = sensitivity (averaged over all spike types)...
h(1) = subplot(1,2,1); % **************
ylabel('Fraction of Actual Spikes'); xlabel('\eta') %%%%%%%%%%%%%%%%%%%**************
patch(xo+w*[0 1 1 0],yo+m_sens*[0 0 1 1],'green');
patch(xo+w*[0 1 1 0],yo+m_sens+m_wr*[0 0 1 1],'red');             % f_wr
patch(xo+w*[0 1 1 0],yo+m_sens+m_wr+m_miss*[0 0 1 1],'white');    % f_miss
if tswitch == 1,
    text(xo+w/8, -0.05, 'actual spikes');
    if m_sens>0.05, text(xo+w/3.2,yo+m_sens/2,'correct (sensitivity)'); end
    if m_wr>0.05, text(xo+w/3,yo+m_sens+m_wr/2,'wrong'); end
    if m_miss>0.05, text(xo+w/5.5,yo+1-m_miss/2,'miss (type II)'); end
end
legend('Correct','Wrong','Missed (Type II)','Location','southwest'); % *********
% text(xo+w/8, -0.05, 'Algorithm'); %%%%%%%%%%%%%%%%%%
set(gca,'fontsize',18)

if N>1
  s_sens = std(sens)/sqrt(N)*sqrt(N - 1);  % unbiased estimate of population std dev
  s_miss = std(f_miss)/sqrt(N)*sqrt(N - 1);
  hold on;
  ylim([0 1])
%   errorbar(xo+w/2*[1 1], yo+[m_sens m_sens+m_wr], [s_sens s_miss],'k.','markersize',10,'linewidth',2);
  terrorbar(xo+w/2*[1 1], yo+[m_sens m_sens+m_wr], [s_sens s_miss], [s_sens s_miss], w*.75, 'units');
end
% set(gca,'xtick',[]); %%%%%%%%%%%%%%%%%%%%%%%%%

% second bar...
h(2) = subplot(1,2,2); %******************
ylabel('Fraction of Detected Spikes'); xlabel('\eta'); %%%%%%%%%%%%%%%%%%*****************
% xo = xo+0.5;  % shift to right %************************
patch(xo+w*[0 1 1 0],yo+m_cor*[0 0 1 1],'green');                 % f_cor
patch(xo+w*[0 1 1 0],yo+m_cor+m_wrd*[0 0 1 1],'red');             % f_wrd
patch(xo+w*[0 1 1 0],yo+m_cor+m_wrd+m_fp*[0 0 1 1],[0 .9 1]);     % f_fp
if tswitch == 1,
    text(xo+w/20, -0.05, 'detected spikes');
    if m_cor>0.05, text(xo+w/3.3,yo+m_cor/2,'correct (precision)'); end
    if m_wrd>0.05, text(xo+w/3.1,yo+m_cor+m_wrd/2,'wrong'); end
    if m_fp>0.05, text(xo+w/16,yo+1-m_fp/2,'false positive (type I)'); end
end
legend('Correct','Wrong','False Pos. (Type I)','Location','southwest'); % ***********
% text(xo+w/8, -0.05, 'Algorithm'); %%%%%%%%%%%%%%%

if N>1
  s_cor = std(f_cor)/sqrt(N)*sqrt(N - 1);  % unbiased estimate of population std dev
  s_fp = std(f_fp)/sqrt(N)*sqrt(N - 1);
  hold on;
  ylim([0 1])
%   errorbar(xo+w/2*[1 1], yo+[m_cor m_cor+m_wrd], [s_cor s_fp],'k.','markersize',10,'linewidth',2);
  terrorbar(xo+w/2*[1 1], yo+[m_cor m_cor+m_wrd], [s_cor s_fp], [s_cor s_fp], w*.75, 'units');
end

% set(gca,'xtick',[]); %%%%%%%%%%%%%%%%%%%%%%%%*********************
set(gca,'fontsize',18)
linkaxes(h,'xy'); %*************

%text(0,-0.1,sprintf('sensitivity = %2d%%\n',100*sens))
% specificity can't be known from this data for each type
% since would need to know what type false positive each was.
%%%%%


function demo_paulabars
% figure;
% paulabars(0,0,10,7,80,100,20,1); title('no errorbars');
% axis equal
a=20;  % simulate some variation of this size
n=10;  % # run batches
figure;
paulabars(2,0,10+randi(a,1,n),7+randi(a,1,n),80+randi(a,1,n),100+randi(a,1,n),20+randi(a,1,n),1);
% title('with errorbars');
axis equal