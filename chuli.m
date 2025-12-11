

%% 0.!!!需手动调整参数
data_path = 'F:\python\EEG\data\d2\';
svae_path = 'F:\python\EEG\data\d2\';
file_name_target = 'sub_';
file_name_nontarget = 'Normal';
stuct_name =  'Data';

raw_temp_data = load ([data_path ,file_name_target,'001']);
raw_temp_data = raw_temp_data.data;
fs_raw = 500;
fs_down = 256; %降采样后采样率
subject_num = [1;122];
temp_abandon = 20; %去掉的暂态，以秒为单位
epoch_length = 512;
perocess_channel = [1 ; 64];
%% 1.Standard_input
% 0.1 rest_data 输入数据格式为：通道数*采样点数
standard_time = floor(size(raw_temp_data,2)/fs_raw- temp_abandon);
remain_raw_time = standard_time - temp_abandon;
remain_trial = floor(remain_raw_time/epoch_length);
remain_confrim_time = remain_trial * epoch_length;

disp(['数据总时长: ' , num2str(standard_time) , '||暂态时长: ' , num2str(temp_abandon), '||修正剩余时长: ' , num2str(remain_confrim_time)]);
disp(['采样率: ' , num2str(fs_raw), '||每段时长: ' , num2str(epoch_length), '||试次数量: ' , num2str(remain_trial)]);
%% 1.切分数据
disp(['数据分段中...']);
Standard_input_target = [];
Standard_input_nontarget = [];

Standard_input_target.fs = fs_down;
Standard_input_nontarget.fs = fs_down;

Standard_input_target.subject_num = subject_num;
Standard_input_nontarget.subject_num = subject_num;

[Standard_input_target.data,Standard_input_nontarget.data] = Standard_cut_rest(data_path,stuct_name,file_name_target,file_name_nontarget,fs_raw,fs_down,subject_num,remain_trial,epoch_length,temp_abandon,perocess_channel);
Standard_input_target.file = file_name_target;
Standard_input_nontarget.file = file_name_nontarget;
disp(['||已完成标准分段||']);

%% 2.保存标准输入文件
disp(['标准分段保存中...']);
save([ svae_path , 'Standard_input_target_',num2str(subject_num(1,1)),'_',num2str(subject_num(2,1))],'Standard_input_target');
save([ svae_path , 'Standard_input_nontarget_',num2str(subject_num(1,1)),'_',num2str(subject_num(2,1))],'Standard_input_nontarget');
disp(['||已完成标准分段保存||']);
t_Standard_input_cost = toc;
disp(['标准输入格式调整完毕，耗时： ',num2str(t_Standard_input_cost)]);
clc;
clear;
function[Standard_input_target,Standard_input_nontarget] = Standard_cut_rest(data_path,stuct_name,file_name_target,file_name_nontarget,fs_raw,fs_down,subject_num,remain_trial,epoch_length,temp_abandon,perocess_channel)
%% 0.参数说明
% data_path 数据路径
% file_name_target 目标数据统一名
% file_name_nontarget 非目标数据统一名
% fs_raw 原始数据采样率
% subject_num 被试起止编号
% remain_trial 剩余试次数
% epoch_length 单试次时长
% temp_abandon 暂态去除时长
% perocess_channel 有效数据通道的首末
Standard_input_target = cell(remain_trial,subject_num(2,1)-subject_num(1,1)+1);
Standard_input_nontarget = cell(remain_trial,subject_num(2,1)-subject_num(1,1)+1);

%% 1.切分数据
a=stuct_name;
b=file_name_nontarget;
print_count=0;
sub_count = 1;
for sub_loop = subject_num(1,1):subject_num(2,1)

    raw_data_target = load ([data_path ,file_name_target , sprintf('%03d', sub_loop)]);
    raw_data_target = raw_data_target.data;
    raw_data_target = downsample(raw_data_target',fs_raw/fs_down)';

    standard_time = floor(size(raw_data_target,2)/fs_raw- temp_abandon);
    remain_raw_time = standard_time - temp_abandon;
    remain_trial = floor(remain_raw_time/epoch_length);
    remain_confrim_time = remain_trial * epoch_length;


    temp_abandon_point = temp_abandon * fs_down;
    for trial_loop = 1:remain_trial
        cut_temp = [];
        cut_temp = raw_data_target(perocess_channel(1,1):perocess_channel(2,1), temp_abandon_point+ (trial_loop-1)*fs_down*epoch_length+1 : temp_abandon_point + trial_loop*fs_down*epoch_length);
        Standard_input_target{trial_loop,sub_count} = cut_temp;

    end
    % disp(['输入标准化分段: ' , num2str(sub_loop/subject_num)]);

    fprintf(repmat('\b',1,print_count));   
    print_count=fprintf('计算标准化分段进度 : %f',sub_count/(subject_num(2,1) - subject_num(1,1)+1));
    sub_count = sub_count+1;
end
fprintf('\n');
end



