
%% 0.特征候选集-参数设置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_path = 'F:\python\EEG\data\ds003474-download\sub-001\eeg\';
svae_path = 'F:\python\EEG\data\ds003474-download\sub-001\eeg\';

subject_num = [1;122];

% Featute_content = ['time\','freq\','time_freq\','space'];
Featute_content = ['time\','freq\','time_freq\'];
% Featute_time_content = ['cross_zero\','std\','apen\','sampentropy\','ar\'];
Featute_time_content = ['cross_zero\','std\','ar\'];
disp(['||特征候选集-参数设置||']);
disp(['特征域内容: ' , Featute_content]);
disp(['时域-候选集: ' , Featute_time_content]);


%% 1.标准输入赋值
Proprocess_target_file = load([data_path ,'Proprocess_target_',num2str(subject_num(1,1)),'_',num2str(subject_num(2,1))]);

stuct_target_name =  'Proprocess_target';
stuct_nontarget_name =  'Proprocess_nontarget';

Proprocess_target_data = Proprocess_target_file.(stuct_target_name).data;


subject_num = Proprocess_target_file.(stuct_target_name).subject_num;
fs_down = Proprocess_target_file.(stuct_target_name).fs_down;

%remain_trial_target = Proprocess_target_file.(stuct_target_name).remain_trial;
remain_trial_target=Proprocess_trial_remain(Proprocess_target_data);

disp(['目标试次剩余: ' , num2str(remain_trial_target),'||平均： ', num2str(mean(remain_trial_target))]);


%% 2.时域特征候选集
if contains(Featute_content,'time')
disp(['时域特征计算中...']);
tic;
[Festure_time_target,Festure_time_candidate_num_target]= Festure_candidate_time(Proprocess_target_data,Featute_time_content,remain_trial_target);    
  
t_time_candidate_cost = toc;

Festure_candidate_time_target = [];
Festure_candidate_time_target.data  = Festure_time_target;
Festure_candidate_time_target.Featute_time_content  = Featute_time_content;
Festure_candidate_time_target.remain_trial_target  = remain_trial_target;
Festure_candidate_time_target.Festure_time_candidate_num_target  = Festure_time_candidate_num_target;
Festure_candidate_time_target.fs_down = fs_down;


disp(['时域特征计算完毕，耗时： ',num2str(t_time_candidate_cost)]);
disp(['时域特征保存中...']);
save([ svae_path , 'Festure_candidate_time_target_',num2str(subject_num(1,1)),'_',num2str(subject_num(2,1))],'Festure_candidate_time_target');
disp(['时域特征保存完毕']);
end


function saen = SampEntropy( dim, r, data, tau )
% SAMPEN Sample Entropy
%   calculates the sample entropy of a given time series data

%   SampEn is conceptually similar to approximate entropy (ApEn), but has
%   following differences:
%       1) SampEn does not count self-matching. The possible trouble of
%       having log(0) is avoided by taking logarithm at the latest step.
%       2) SampEn does not depend on the datasize as much as ApEn does. The
%       comparison is shown in the graph that is uploaded.

%   dim     : embedded dimension
%   r       : tolerance (typically 0.2 * std)
%   data    : time-series data
%   tau     : delay time for downsampling (user can omit this, in which case
%             the default value is 1)
%

if nargin < 4, tau = 1; end
if tau > 1, data = downsample(data, tau); end

N = length(data);
result = zeros(1,2);

for m = dim:dim+1
    Bi = zeros(1,N-m+1);
    dataMat = zeros(m,N-m+1);

    % setting up data matrix
    for i = 1:m
        dataMat(i,:) = data(i:N-m+i);
    end

    % counting similar patterns using distance calculation
    for j = 1:N-m+1
        % calculate Chebyshev distance, excluding self-matching case
        dist = max(abs(dataMat - repmat(dataMat(:,j),1,N-m+1)));
        % calculate Heaviside function of the distance
        % User can change it to any other function
        % for modified sample entropy (mSampEn) calculation
        D = (dist <= r);
        % excluding self-matching case
        Bi(j) = (sum(D)-1)/(N-m);
    end

    % summing over the counts
    result(m-dim+1) = sum(Bi)/(N-m+1);

end

saen = -log(result(2)/result(1));

end




function apen = ApEn( dim, r, data, tau )
%ApEn
%   dim : embedded dimension
%   r : tolerance (typically 0.2 * std)
%   data : time-series data
%   tau : delay time for downsampling

%   Changes in version 1
%       Ver 0 had a minor error in the final step of calculating ApEn
%       because it took logarithm after summation of phi's.
%       In Ver 1, I restored the definition according to original paper's
%       definition, to be consistent with most of the work in the
%       literature. Note that this definition won't work for Sample
%       Entropy which doesn't count self-matching case, because the count
%       can be zero and logarithm can fail.
%
%       A new parameter tau is added in the input argument list, so the users
%       can apply ApEn on downsampled data by skipping by tau.
%---------------------------------------------------------------------
% coded by Kijoon Lee,  kjlee@ntu.edu.sg
% Ver 0 : Aug 4th, 2011
% Ver 1 : Mar 21st, 2012
%---------------------------------------------------------------------
if nargin < 4, tau = 1; end
if tau > 1, data = downsample(data, tau); end

N = length(data);
result = zeros(1,2);

for j = 1:2
    m = dim+j-1;
    phi = zeros(1,N-m+1);
    dataMat = zeros(m,N-m+1);

    % setting up data matrix
    for i = 1:m
        dataMat(i,:) = data(i:N-m+i);
    end

    % counting similar patterns using distance calculation
    for i = 1:N-m+1
        tempMat = abs(dataMat - repmat(dataMat(:,i),1,N-m+1));
        boolMat = any( (tempMat > r),1);
        phi(i) = sum(~boolMat)/(N-m+1);
    end

    % summing over the counts
    result(j) = sum(log(phi))/(N-m+1);
end

apen = result(1)-result(2);

end



function [Festure_time,Festure_time_candidate_num]= Festure_candidate_time(Standard_input_data,Featute_time_content,remain_trial)

Festure_time = [];
%% 1.cross_zero
cross_zero = [];
if contains(Featute_time_content,'cross_zero')
cross_zero = zeros(sum(remain_trial),size(Standard_input_data{1,1},1));
count_trial = 1;
for sub_loop = 1:size(remain_trial,2)
for trial_loop = 1:remain_trial(1,sub_loop)
cross_zero_channel_temp = zeros(1,size(Standard_input_data{1,1},1));
for channel_loop = 1:size(Standard_input_data{1,1},1)
    for point_loop = 1:size(Standard_input_data{1,1},2)-1
    if Standard_input_data{trial_loop,sub_loop}(channel_loop,point_loop) * Standard_input_data{trial_loop,sub_loop}(channel_loop,point_loop+1)<0
    cross_zero_channel_temp(1,channel_loop) = cross_zero_channel_temp(1,channel_loop) + 1;
    end
    end
end
cross_zero(count_trial,:) = cross_zero_channel_temp;   
count_trial = count_trial+1;
end
end
end

%% 2.std
fest_std = [];
if contains(Featute_time_content,'std')
fest_std = zeros(sum(remain_trial),size(Standard_input_data{1,1},1));
count_trial = 1;
for sub_loop = 1:size(remain_trial,2)
for trial_loop = 1:remain_trial(1,sub_loop)
std_temp  = [];
std_temp  = std(Standard_input_data{trial_loop,sub_loop}');
fest_std(count_trial,:) = std_temp;
count_trial = count_trial+1;    
end
end
end

%% 3.近似熵 
fest_apen = [];
if contains(Featute_time_content,'apen')
r_apen = 0.2*fest_std';    
fest_apen_2 = zeros(sum(remain_trial),size(Standard_input_data{1,1},1));
fest_apen_3 = zeros(sum(remain_trial),size(Standard_input_data{1,1},1));
count_trial = 1;

for sub_loop = 1:size(remain_trial,2)
for trial_loop = 1:remain_trial(1,sub_loop)
for channel_loop = 1:size(Standard_input_data{1,1},1)
fest_apen_2(count_trial,channel_loop) = ApEn( 2, r_apen(1,count_trial), Standard_input_data{trial_loop,sub_loop}(channel_loop,:), 1 );
fest_apen_3(count_trial,channel_loop) = ApEn( 3, r_apen(1,count_trial), Standard_input_data{trial_loop,sub_loop}(channel_loop,:), 1 );
end
count_trial = count_trial+1;    
end
end

fest_apen = [fest_apen_2 fest_apen_3];
end

%% 4.样本熵 
fest_sampentropy = [];
if contains(Featute_time_content,'sampentropy')
r_sampentropy = 0.2*fest_std';    
fest_sampentropy_2 = zeros(sum(remain_trial),size(Standard_input_data{1,1},1));
fest_sampentropy_3 = zeros(sum(remain_trial),size(Standard_input_data{1,1},1));
count_trial = 1;

for sub_loop = 1:size(remain_trial,2)
for trial_loop = 1:remain_trial(1,sub_loop)
for channel_loop = 1:size(Standard_input_data{1,1},1)
fest_sampentropy_2(count_trial,channel_loop) = SampEntropy( 2, r_sampentropy(1,count_trial), Standard_input_data{trial_loop,sub_loop}(channel_loop,:), 1 );
fest_sampentropy_3(count_trial,channel_loop) = SampEntropy( 3, r_sampentropy(1,count_trial), Standard_input_data{trial_loop,sub_loop}(channel_loop,:), 1 );
end
count_trial = count_trial+1;    
end
end

fest_sampentropy = [fest_sampentropy_2 fest_sampentropy_3];
end


%% 5.AR
fest_ar = [];
if contains(Featute_time_content,'ar')
ar_order = 8;
fest_ar = zeros(sum(remain_trial),size(Standard_input_data{1,1},1)*ar_order);
count_trial = 1;
for sub_loop = 1:size(remain_trial,2)
for trial_loop = 1:remain_trial(1,sub_loop)
ar_temp  = [];
ar_temp = aryule(Standard_input_data{trial_loop,sub_loop}',ar_order);
ar_temp = ar_temp(:,2:ar_order+1);
fest_ar(count_trial,:) = reshape(ar_temp',1,size(Standard_input_data{1,1},1)*ar_order);
count_trial = count_trial+1;    
end
end
end

%% 时域特征合并
Festure_time = [cross_zero   fest_std    fest_apen    fest_sampentropy   fest_ar];
Festure_time_candidate_num = [size(cross_zero,2)   size(fest_std,2)    size(fest_apen,2)    size(fest_sampentropy,2)   size(fest_ar,2)];
end


