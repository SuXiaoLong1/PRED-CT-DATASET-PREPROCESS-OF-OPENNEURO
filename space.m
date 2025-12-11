
%% 0.特征候选集-参数设置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_path = 'F:\python\EEG\data\ds003474-download\sub-001\eeg\';
svae_path = 'F:\python\EEG\data\ds003474-download\sub-001\eeg\';

subject_num = [1;122];

Featute_space_content = ['csp\'];
space_filter_num = 20;
Featute_content = ['space'];



%% 1.标准输入赋值
Proprocess_target_file = load([data_path ,'Proprocess_target_',num2str(subject_num(1,1)),'_',num2str(subject_num(2,1))]);

stuct_target_name =  'Proprocess_target';


Proprocess_target_data = Proprocess_target_file.(stuct_target_name).data;


subject_num = Proprocess_target_file.(stuct_target_name).subject_num;
fs_down = Proprocess_target_file.(stuct_target_name).fs_down;

%remain_trial_target = Proprocess_target_file.(stuct_target_name).remain_trial;
remain_trial_target=Proprocess_trial_remain(Proprocess_target_data);


disp(['目标试次剩余: ' , num2str(remain_trial_target),'||平均： ', num2str(mean(remain_trial_target))]);



%% 5.空域特征候选集
if contains(Featute_content,'space')
disp(['空域特征计算中...']);

tic;
[Festure_space_target,Festure_space_candidate_num_target]= Festure_candidate_space(Proprocess_target_data,remain_trial_target,Featute_space_content,space_filter_num);
t_space_candidate_cost = toc;
disp(['空域特征计算完毕，耗时(秒)： ',num2str(t_space_candidate_cost)]);

Festure_candidate_space_target = [];
Festure_candidate_space_target.data  = Festure_space_target;
Festure_candidate_space_target.Featute_space_content  = Featute_space_content;
Festure_candidate_space_target.remain_trial_target  = remain_trial_target;
Festure_candidate_space_target.Festure_space_candidate_num_target  = Festure_space_candidate_num_target;
Festure_candidate_space_target.fs_down = fs_down;



disp(['空域特征保存中...']);
save([ svae_path , 'Festure_candidate_space_target_',num2str(subject_num(1,1)),'_',num2str(subject_num(2,1))],'Festure_candidate_space_target');
save([ svae_path , 'Festure_candidate_space_nontarget_',num2str(subject_num(1,1)),'_',num2str(subject_num(2,1))],'Festure_candidate_space_nontarget');
disp(['空域特征保存完毕']);

end

function [Pd Pf]=ROC( Y_taris,Y_othis)
amin=min(min(Y_othis));
amax=max(max(Y_othis));
middata=(amax-amin)/1e5;
x=[amin:middata:amax];   %x 表示阈值选择范围
Pd=zeros(1,length(x));
Pf=zeros(1,length(x));
for i=1:length(x)
   Pd(i)=length(find((Y_taris>x(i))==1))/length(Y_taris);
   Pf(i)=length(find((Y_othis>x(i))==1))/length(Y_othis);
end
% figure
% plot(Pf,Pd)
% for i=1:length(x)
%    dt(i)=Pd(i)/Pf(i) ;
% end
%   b=x(find(dt==max(dt)));
end


function value = AUC(Pf,Pd)
% 给定Pf，Pd返回对应的AUC值
% s=size(Pf,2);
% value=sum(Pd)*(1/s);
Pf_diff=diff(Pf);
value=abs(sum(Pd(2:end).*Pf_diff));
end



function [Festure_space_target,Festure_space_candidate_num_target]= Festure_candidate_space(Proprocess_target_data,remain_trial_target,Featute_space_content,space_filter_num)
%% 注：空间特征属于需联合其他特征域的二级特征
%% 注：空间特征一般不处理原始数据，如CSP一般在数据分析后，结合时频的差异性特征，构建时空 or 空频特征

% https://cloud.tencent.com/developer/article/1654056 源代码网址
% CSp的推导 https://blog.csdn.net/MissXy_/article/details/81264953

% 输入参数
% Proprocess_target_data 目标数据（预处理后）
% Proprocess_nontarget_data 非目标数据（预处理后）
% remain_trial_target 目标试次剩余数量
% remain_trial_nontarget 非目标试次剩余数量
% Featute_space_content 空域特征内容
% space_filter_num 一般的CSP特征个数
% 
% 输出参数
% Festure_space_target 目标的空间特征
% Festure_space_candidate_num_target 目标各类的空间特征个数
% Festure_space_nontarget 非目标的非空间特征
% Festure_space_candidate_num_nontarget 非目标各类的空间特征个数

%% 1.数据格式整合
target_3d = [];
trial_count = 1;
for sub_loop = 1:size(Proprocess_target_data,2)
for trial_loop = 1:remain_trial_target(1,sub_loop)
target_3d(:,:,trial_count) = Proprocess_target_data{trial_loop,sub_loop}';
trial_count = trial_count + 1;
end
end



EEG_space_data = target_3d;
EEG_space_label = [zeros(sum(remain_trial_target),1)];

%% 2\csp特征计算
if contains(Featute_space_content,'csp')
%check and initializations
EEG_Channels = size(EEG_space_data,2);
EEG_Trials = size(EEG_space_data,3);
classLabels = unique(EEG_space_label);% Return non-repeating values
EEG_Classes = length(classLabels);
covMatrix = cell(EEG_Classes,1); % 协方差矩阵
% Computing the normalized covariance matrices for each trial
trialCov = zeros(EEG_Channels,EEG_Channels,EEG_Trials);
for i = 1:EEG_Trials
    E = EEG_space_data(:,:,i)';
    EE = E*E';
    trialCov(:,:,i) = EE./trace(EE);  % 计算协方差矩阵
end
clear E;
clear EE;
% 计算每一类样本数据的空间协方差之和
for i = 1:EEG_Classes
    covMatrix{i} = mean(trialCov(:,:,EEG_space_label == classLabels(i)),3);
end
% 计算两类数据的空间协方差之和
covTotal = covMatrix{1} + covMatrix{2};
% 计算特征向量和特征矩阵
[Uc,Dt] = eig(covTotal);
% 特征值要降序排列
eigenvalues = diag(Dt);
[eigenvalues,egIndex] = sort(eigenvalues, 'descend');% 降序
Ut = Uc(:,egIndex);
% 矩阵白化
P = diag(sqrt(1./eigenvalues))*Ut';
% 矩阵P作用求公共特征向量transformedCov1 
transformedCov1 = P*covMatrix{1}*P';
%计算公共特征向量transformedCov1的特征向量和特征矩阵
[U1,D1] = eig(transformedCov1);
eigenvalues = diag(D1);
[eigenvalues,egIndex] = sort(eigenvalues, 'descend');% 降序排列
U1 = U1(:, egIndex);
% 计算投影矩阵W
CSPMatrix = U1' * P;
% 计算特征矩阵
% CSP特征选择参数m    CSP特征为2*m个
features_train = zeros(EEG_Trials, 2*space_filter_num+1);
features_test = zeros(EEG_Trials, 2*space_filter_num+1);
Filter = CSPMatrix([1:space_filter_num (end-space_filter_num+1):end],:);
%extracting the CSP features from each trial
for t=1:EEG_Trials    
    %projecting the data onto the CSP filters    
    projectedTrial_train = Filter * EEG_space_data(:,:,t)';    
    %generating the features as the log variance of the projected signals
    variances_train = var(projectedTrial_train,0,2);  
    for f=1:length(variances_train)
        features_train(t,f) = log(variances_train(f));
        % features_train(t,f) = log(variances_train(f)/sum(variances_train));   %修改后对应公式
    end

end
CSP_Train_feature = features_train(:,1:2*space_filter_num);
end



Festure_space_target  = CSP_Train_feature(1:sum(remain_trial_target),:);
Festure_space_candidate_num_target = space_filter_num*2;


end

function [w,v,dev,stats]=HDCA_train(w_blocknum,train_data,train_label)
ch=size(train_data,1);

 target_id = find(train_label==1);
 other_id = find(train_label==0);
 tar_num=length(target_id);
 oth_num=length(other_id);

 target_signal=train_data(:,:,target_id);
 other_signal=train_data(:,:,other_id);

%  w_blocknum=floor(size(train_data,2)/fs*40);               %在计算w权重时对波形的分段,即认为1s的数据分为w_blocknum段，每段1000/w_blocknum ms
 w=zeros(ch,w_blocknum);
 num=floor(size(train_data,2)/w_blocknum);          %计算每一分段多少点，要求处下来必须为整数
%  Y_tar_sig=zeros(tar_num,200);
%  Y_oth_sig=zeros(oth_num,200);
 %% 对16导数据进行分块处理
tic
 for i=1:w_blocknum                 %计算每个分段权重，并将16导联数据分段加权到一个信号上
    tar_blockms=target_signal(:,num*(i-1)+1:num*i,:);
    tar_blockms2=mean(tar_blockms,2);        %16导块平均部分
    tar_reshape=reshape(tar_blockms2,ch,tar_num);
    oth_blockms=other_signal(:,num*(i-1)+1:num*i,:);
    oth_blockms2=mean(oth_blockms,2);
    oth_reshape=reshape(oth_blockms2,ch,oth_num); 
    [w(:,i),~]=Fisher(tar_reshape,oth_reshape);
    for t=1:tar_num
        Y_tar_sig(t,num*(i-1)+1:num*i)=(w(:,i)')*tar_blockms(:,:,t);
    end
    for t=1:oth_num
        Y_oth_sig(t,num*(i-1)+1:num*i)=(w(:,i)')*oth_blockms(:,:,t);
    end
 end
 
 %% 取分块间平均值
  v_blocknum=w_blocknum;
 Y_tar_blomean=zeros(tar_num,v_blocknum);
 Y_oth_blomean=zeros(oth_num,v_blocknum);

 for i=1:tar_num
     Y_tar_blomean(i,:)=block_mean(Y_tar_sig(i,:),v_blocknum);
 end
 for i=1:oth_num
    Y_oth_blomean(i,:)= block_mean(Y_oth_sig(i,:),v_blocknum);
 end
%%
% train_data2=cat(1,Y_tar_sig,Y_oth_sig);
% pca_rawdata=reshape(train_data2,tar_num+oth_num,300);
%   [coeff,score1,latent] = pca(pca_rawdata);
%    PcaData=score1(:,1:50)';
%  data=reshape(PcaData,30,tar_num+oth_num);
%%
%  save('BPtrain.mat','Y_tar_blomean','Y_oth_blomean','Y_tar_sig','Y_oth_sig');
 %% 训练时间权重v
% %  [v,b]=Fisher(Y_tar_blomean',Y_oth_blomean');
 train_data=[Y_tar_blomean;Y_oth_blomean];
%  train_data=PcaData;
 train_label=[ones(size(Y_tar_sig,1),1);zeros(size(Y_oth_sig,1),1)];
 [v,dev,stats]=glmfit(train_data,train_label,'binomial','link','logit');

end

