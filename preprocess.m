
%% 0.标准数据参数设置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_path = 'F:\python\EEG\data\d2\';
svae_path = 'F:\python\EEG\data\d2\';

proprocess_content = ['baseline\','filter\','channel_abandon\','auto_abandon\','trial_abandon'];
channel_abandon_num = [];

trial_threshold = 100;
subject_num = [1 ; 122];
Baseline_reference = [0;0.2];
filter_low_para = [80;100];
filter_high_para = [0.01;1];

outlier_threshold = 3;

disp(['||预处理参数设置||']);
disp(['预处理内容: ' , proprocess_content]);
disp(['基线长度: ' , num2str(Baseline_reference(1,1)),'-',num2str(Baseline_reference(2,1))]);
disp(['低通起止: ' , num2str(filter_low_para(1,1)),'-',num2str(filter_low_para(2,1))]);
disp(['高通起止: ' , num2str(filter_high_para(1,1)),'-',num2str(filter_high_para(2,1))]);
disp(['通道舍弃: ' , num2str(channel_abandon_num')]);
disp(['试次幅度阈值: ' , num2str(trial_threshold)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.标准输入赋值
disp(['||1.静息态or任务态-标准输入数据加载中...||']);

Standard_target_file = load([data_path , 'Standard_input_target_',num2str(subject_num(1,1)),'_',num2str(subject_num(2,1))]);
stuct_target_name =  'Standard_input_target';



Standard_target_data = Standard_target_file.(stuct_target_name).data;

subject_num = Standard_target_file.(stuct_target_name).subject_num;
fs_down = Standard_target_file.(stuct_target_name).fs;
trial_every_sub = size(Standard_target_data,1);
disp(['被试量: ' , num2str(subject_num(1,1)),'-',num2str(subject_num(2,1))]);

if (filter_low_para(2,1)>fs_down/2)
    disp(['低通滤波参数不符合奈奎斯特带宽，请调高fs_down或降低低通参数']);
end
%% 2.预处理
%2.1 基线矫正
if contains(proprocess_content,'baseline')
Proprocess_baseline_target = Proprocess_baseline(Standard_target_data,fs_down,Baseline_reference);

end
%2.2 滤波
if contains(proprocess_content,'filter')
Proprocess_filter_target = Proprocess_filter(Proprocess_baseline_target,fs_down,filter_low_para,filter_high_para);

end
%2.3 多余通道剔除
if contains(proprocess_content,'channel_abandon')
Proprocess_channel_abandon_target = Proprocess_channel_abandon(Proprocess_filter_target,channel_abandon_num);

end
%2.4 异常通道自动替换
if contains(proprocess_content,'auto_abandon')
[Proprocess_auto_abandon_target,auto_channel_list_target] = Proprocess_auto_abandon(Proprocess_channel_abandon_target,outlier_threshold);

end
%2.5 试次剔除
if contains(proprocess_content,'trial_abandon')
Proprocess_trial_abandon_target = Proprocess_trial_abandon(Proprocess_auto_abandon_target,trial_threshold);

end


%[remain_target_trial,remain_nontarget_trial]= Proprocess_trial_remain(Proprocess_trial_abandon_target,Proprocess_trial_abandon_nontarget);
% disp(['目标试次剩余比例: ' , num2str(remain_target_trial/size(Proprocess_trial_abandon_target,1)),'||平均： ', num2str(mean(remain_target_trial/size(Proprocess_trial_abandon_target,1)))]);
% disp(['非目标试次剩余比例: ' , num2str(remain_nontarget_trial/size(Proprocess_trial_abandon_nontarget,1)),'||平均： ', num2str(mean(remain_nontarget_trial/size(Proprocess_trial_abandon_nontarget,1)))]);
remain_target_trial = Proprocess_trial_remain(Proprocess_trial_abandon_target);
disp(['目标试次剩余比例: ' , num2str(remain_target_trial/size(Proprocess_trial_abandon_target,1)),'||平均： ', num2str(mean(remain_target_trial/size(Proprocess_trial_abandon_target,1)))]);
%% 3.预处理数据保存
Proprocess_target = [];
Proprocess_target.remain_trial = remain_target_trial;
Proprocess_target.fs_down = fs_down;
Proprocess_target.subject_num = subject_num;
Proprocess_target.data = Proprocess_trial_abandon_target;
Proprocess_target.Baseline_reference = Baseline_reference;


disp(['标准分段保存中...']);
save([ svae_path , 'Proprocess_target_',num2str(subject_num(1,1)),'_',num2str(subject_num(2,1))],'Proprocess_target');
disp(['||已完成标准分段保存||']);


function baseline_out= Proprocess_baseline(Standard_input_data,fs_down,Baseline_reference)
% Standard_input_data 标准输入的cell脑电数据,cell(试次数*被试数)[通道数，时间点数]
% fs_down 标准输入时降采样率
% Baseline_reference 极限矫正的参考时间段，一般为0~0.2的试次前均值，或者0~1的全时段均值

baseline_start = floor(Baseline_reference(1,1)*fs_down)+1;
baseline_end = floor(Baseline_reference(2,1)*fs_down);

baseline_out = cell(size(Standard_input_data,1),size(Standard_input_data,2));

for trial_loop = 1:size(Standard_input_data,1)
for sub_loop = 1:size(Standard_input_data,2)
for channel_loop = 1:size(Standard_input_data{1,1},1)
    if ~isempty(Standard_input_data{trial_loop,sub_loop})
    baseline_out{trial_loop,sub_loop}(channel_loop,:) = Standard_input_data{trial_loop,sub_loop}(channel_loop,:) - mean(Standard_input_data{trial_loop,sub_loop}(channel_loop,baseline_start:baseline_end));    
    end
end
end
end

end


function [filter_out]= Proprocess_filter(Standard_input_data,fs_down,filter_low_para,filter_high_para)

%导入参数
low_start = filter_low_para(1,1);
low_end = filter_low_para(2,1);
high_start =filter_high_para(1,1);
high_end = filter_high_para(2,1);
filter_out = cell(size(Standard_input_data,1),size(Standard_input_data,2));
% 低通
Rp_low = 0.5;
Rs_low =5;
[N_low,Wpo_low]=cheb1ord(2*low_start/fs_down,2*low_end/fs_down,Rp_low,Rs_low);
[b_low,a_low]=cheby1(N_low,Rp_low,Wpo_low,'low');

% 高通
Rp_high =1;
Rs_high =10;
[N_high,Wpo_high]=cheb1ord(2*high_end/fs_down,2*high_start/fs_down,Rp_high,Rs_high);
[b_high,a_high]=cheby1(N_high,Rp_high,Wpo_high,'high');

% 滤波
for trial_loop = 1:size(Standard_input_data,1)
for sub_loop = 1:size(Standard_input_data,2)
for channel_loop = 1:size(Standard_input_data{1,1},1)    
    if ~isempty(Standard_input_data{trial_loop,sub_loop})
    
    temp_filter = [];
    temp_filter = filter(b_low,a_low,Standard_input_data{trial_loop,sub_loop}(channel_loop,:));
    temp_filter = filter(b_high,a_high,temp_filter);     
    filter_out{trial_loop,sub_loop}(channel_loop,:) = temp_filter;  
    end
end
end
end

end


function channel_abandon_out= Proprocess_channel_abandon(Standard_input_data,channel_abandon_num)
% Standard_input_data 标准输入的cell脑电数据,cell(试次数*被试数)[通道数，时间点数]
% channel_abandon_num 舍弃通道编号
abandon_size = size(channel_abandon_num,1);

channel_abandon_out = cell(size(Standard_input_data,1),size(Standard_input_data,2));

for sub_loop = 1:size(Standard_input_data,2)
for trial_loop = 1:size(Standard_input_data,1)
    channel_count = 1;
for channel_loop = 1:size(Standard_input_data{1,1},1)
if ~isempty(Standard_input_data{trial_loop,sub_loop})  
if isempty(find(channel_abandon_num==channel_loop))
    channel_abandon_out{trial_loop,sub_loop}(channel_count,:) = Standard_input_data{trial_loop,sub_loop}(channel_loop,:);    
    channel_count = channel_count+1;
end
end
end
end
end

end


function [auto_abandon_data,auto_channel_list]= Proprocess_auto_abandon(Standard_input_data,outlier_threshold)
%% 注：auto_abandon只是将抖动异常通道使用相邻通道替换，而不影响数据通道个数，处理后数据通道个数仍一致

% Standard_input_data 标准输入的cell脑电数据,cell(试次数*被试数)[通道数，时间点数]
% outlier_threshold 为异常点的阈值，即异常点抖动是均值的几倍则被判为异常点，一般设为5
% auto_abandon_data 为自动修正（用周围通道替换）通道后的数据
% auto_channel_list 为各被试替换的通道列表

channel_temp_std = [];
auto_channel_list = cell(1,size(Standard_input_data,2));
for sub_loop = 1:size(Standard_input_data,2)
    abandon_count = 1;
    channel_temp_std = [];
    for trial_loop = 1:size(Standard_input_data,1)
        channel_temp_std(:,trial_loop) = std(Standard_input_data{trial_loop,sub_loop}')';
    end   
    std_mean = mean(channel_temp_std');
    abandon_level = mean(outlier_threshold*std_mean);
    for channel_loop = 1:size(Standard_input_data{1,1},1)
        if std_mean(1,channel_loop) > abandon_level
            auto_channel_list{1,sub_loop}(1,abandon_count) = channel_loop;
        end
    end
end
    
    
for sub_loop = 1:size(Standard_input_data,2)
    if ~isempty(auto_channel_list{1,sub_loop})
    for trial_loop = 1:size(Standard_input_data,1)
        for channel_loop = 1:size(Standard_input_data{1,1},1)
        if ismember(auto_channel_list{1,sub_loop}(1,abandon_count),channel_loop)
            replace_data=[];
            if channel_loop==1
                replace_data = Standard_input_data{trial_loop,sub_loop}(2,:);
            elseif channel_loop==size(Standard_input_data{1,1},1)
                replace_data = Standard_input_data{trial_loop,sub_loop}(size(Standard_input_data{1,1},1)-1,:);
            else
                replace_data = (Standard_input_data{trial_loop,sub_loop}(channel_loop-1,:) + Standard_input_data{trial_loop,sub_loop}(channel_loop+1,:))/2;                
            end
            Standard_input_data{trial_loop,sub_loop}(channel_loop,:) = replace_data;
        end
        end
        end
    end
end

auto_abandon_data = Standard_input_data;

end


function trial_abandon_out= Proprocess_trial_abandon(Standard_input_data,trial_threshold)
% Standard_input_data 标准输入的cell脑电数据,cell(试次数*被试数)[通道数，时间点数]
% trial_abandon_num 舍弃通道编号


trial_abandon_out = cell(size(Standard_input_data,1),size(Standard_input_data,2));

for sub_loop = 1:size(Standard_input_data,2)
    trial_count = 1;
for trial_loop = 1:size(Standard_input_data,1)
if max(max(Standard_input_data{trial_loop,sub_loop}))<trial_threshold
    trial_abandon_out{trial_count,sub_loop}= Standard_input_data{trial_loop,sub_loop};    
    trial_count = trial_count+1;
end
end
end

end

function remain_trial=Proprocess_trial_remain(Proprocess_trial_abandon_target)
%获取每个通道剩余数据量
%分类
remain_trial = zeros(1,size(Proprocess_trial_abandon_target,2));
for i=1:size(Proprocess_trial_abandon_target,2)
   for j=1:size(Proprocess_trial_abandon_target,1)
   if isempty(Proprocess_trial_abandon_target{j,i})==1
        remain_trial(i)=remain_trial(i)+1;
        
   end
   end
end
end
