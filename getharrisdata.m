% make 1-ch of Harris data into matlab format. Barnett 5/19/17
% adapted from /home/alex/SCDA/mountainlab_devel/examples/grab_harris2000_dataset.m 

clear
fname = '/home/alex/SCDA/datasets/Harris2000/d5331/d533101.dat';
fid = fopen(fname,'r');
if fid==-1, error('Harris data not found!'); end
Y = fread(fid,inf,'int16');
fclose(fid);
n = numel(Y);
Nch = 8;
N = n/Nch; if N~=round(N), error('non-integer number of time pts for given M'); 
end
Y = reshape(Y,[Nch N]);
samples_per_sec = 1e4;
t = (1:N)/samples_per_sec;   % time indices
%j = find((t>26 & t<109.49) | t>110.29);   % ahb original cutout excluded opening
j = find(t<109.49 | t>110.29);   % cut out 0.8 sec of bad bursting/noise part, forgets absolute t. As in Ekanadham et al 2013.
N = numel(j);
Y = Y(:,j);

YIC = Y(6,:); % pull out the IC (intra-cellular), for ground truth
trig = (max(YIC)+min(YIC))/2;   % trigger level (checked by eye)
truetimes = 1+find(diff(YIC>trig)==1); % upwards-going times, sample units
%figure; plot(1:N, YIC); vline(d.truetimes); xlabel('t (in samples)');
ntrue =numel(truetimes)

Y = Y(2:5,:);  % the EC channels
Y = Y(1,:); % just the first

addpath /home/alex/SCDA/mountainlab/matlab/processing
Y = ms_bandpass_filter(Y,struct('samplerate',samples_per_sec,'freq_min',300,'freq_max',6000));
Y =Y(1:2392200);  % kill filter spike at end
Y = ms_normalize_channels(Y);
figure; plot(Y);

clipsize=30;
firingtimes = ms_detect3(Y,struct('detect_threshold',5,'detect_interval',10,'clip_size',clipsize,'sign',-1));
nclips = numel(firingtimes)

clips = ms_extract_clips2(Y,firingtimes,clipsize);
[FF, subspace] = ms_event_features(clips,10); % num fea

addpath /home/alex/SCDA/mountainlab/mountainsort/src/isosplit/
addpath /home/alex/SCDA/mountainlab/matlab/msutils/
labels = isosplit2(FF);

K=max(labels)   % hope it's 3, sometimes 2
clips = squeeze(clips);   % 1 ch only
templates = nan(K,clipsize);
for k=1:K, templates(k,:) = mean(clips(:,labels==k),2)'; nk(k)=sum(labels==k); end
nk
figure; plot(templates'); title('mean templates');
figure; for k=1:K, subplot(1,K,k); plot(clips(:,labels==k)); end
Y = single(Y); % save space
sortedlabels=labels;
groundtruthlabels = nan*labels;
for i=1:nclips, if sum(abs(truetimes-firingtimes(i))<5), groundtruthlabels(i)=1; end; end  % crude O(N^2)
groundtruthlabels = groundtruthlabels*2;   % make spike number match
sum(~isnan(groundtruthlabels))  % how many found
save harrisdata.mat K Y clips clipsize firingtimes sortedlabels nclips samples_per_sec templates truetimes groundtruthlabels ntrue
