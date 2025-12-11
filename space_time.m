
%% 0.特征候选集-参数设置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_path = 'F:\python\EEG\data\ds003474-download\sub-001\eeg\';
svae_path = 'F:\python\EEG\data\ds003474-download\sub-001\eeg\';

subject_num = [1;122];
freq_resolut = 2;
freq_scale = [1;60];

Featute_content = ['time\','freq\','time_freq\','space'];

Featute_time_freq_content = ['emd\','DE\'];
imf_level = 2;

disp(['||特征候选集-参数设置||']);
disp(['特征域内容: ' , Featute_content]);



%% 1.标准输入赋值
Proprocess_target_file = load([data_path ,'Proprocess_target_',num2str(subject_num(1,1)),'_',num2str(subject_num(2,1))]);

stuct_target_name =  'Proprocess_target';


Proprocess_target_data = Proprocess_target_file.(stuct_target_name).data;


subject_num = Proprocess_target_file.(stuct_target_name).subject_num;
fs_down = Proprocess_target_file.(stuct_target_name).fs_down;

%remain_trial_target = Proprocess_target_file.(stuct_target_name).remain_trial;
remain_trial_target=Proprocess_trial_remain(Proprocess_target_data);

disp(['目标试次剩余: ' , num2str(remain_trial_target),'||平均： ', num2str(mean(remain_trial_target))]);



%% 4.时频域特征候选集
if contains(Featute_content,'time_freq')
disp(['时频域特征计算中...']);

tic;
[Festure_time_freq_target,Festure_time_freq_candidate_num_target]= Festure_candidate_time_freq(Proprocess_target_data,Featute_time_freq_content,remain_trial_target,fs_down,imf_level);


if contains(Featute_time_freq_content,'DE')
Festure_time_freq_target = log(Festure_time_freq_target);

end
t_time_freq_candidate_cost = toc;
disp(['时频域特征计算完毕，耗时(秒)： ',num2str(t_time_freq_candidate_cost)]);

Festure_candidate_time_freq_target = [];
Festure_candidate_time_freq_target.data  = Festure_time_freq_target;
Festure_candidate_time_freq_target.Featute_time_freq_content  = Featute_time_freq_content;
Festure_candidate_time_freq_target.remain_trial_target  = remain_trial_target;
Festure_candidate_time_freq_target.Festure_time_freq_candidate_num_target  = Festure_time_freq_candidate_num_target;
Festure_candidate_time_freq_target.fs_down = fs_down;



disp(['时频域特征保存中...']);
save([ svae_path , 'Festure_candidate_time_freq_target_',num2str(subject_num(1,1)),'_',num2str(subject_num(2,1))],'Festure_candidate_time_freq_target');

disp(['时频域特征保存完毕']);
end

function five_band_sum = sum_five_band(fft_temp,fs_down)
%% 这只是一行的5频带求和，请在外面加Channe_loop循环

delta =[1;4]; %δ 
theta =[4;8];  %θ
alpha =[8;12];  %α?
beta = [12;30];  %β ? 
gamma =[30;60];  %γ ?
five_band = [delta theta alpha beta gamma];

fft_resolut = fs_down/size(fft_temp,2);
epoch_num = size(five_band,2);
epoch_length = five_band.*fft_resolut;

five_band_sum = [];
for cut_loop = 1:epoch_num
fft_sum_temp = sum(fft_temp(:, epoch_length(1,cut_loop): epoch_length(2,cut_loop)));
five_band_sum = [five_band_sum fft_sum_temp'];
end

end



function [Festure_time_freq,Festure_time_freq_candidate_num]= Festure_candidate_time_freq(Standard_input_data,Featute_time_freq_content,remain_trial,fs_down,imf_level)

Festure_time_freq = [];
%% 1.emd只进行传统5频带 five_band
if contains(Featute_time_freq_content,'emd')
fest_five_band = [];
count_trial = 1;
for sub_loop = 1:size(remain_trial,2)
    disp([sub_loop])
for trial_loop = 1:remain_trial(1,sub_loop)
        b1=[];
        five_band_temp = [];
        for channel_loop = 1:size(Standard_input_data{1,1},1)
            fft_temp = [];
            imf_temp = [];
            imf_temp = emd(Standard_input_data{trial_loop, sub_loop}(channel_loop,:),'Interpolation','pchip','MaxNumIMF',imf_level ,'Display',0);
            if isempty(imf_temp)
                imf_temp=b1;
            end
            fft_temp = abs(fft(imf_temp(:,imf_level)',fs_down));
            five_band_temp(channel_loop,:) =sum_five_band(fft_temp,fs_down);
            b1=imf_temp;
        end

    fest_five_band(count_trial,:) = reshape(five_band_temp',1,size(five_band_temp,1)*size(five_band_temp,2));
    count_trial = count_trial+1;    
end
end
end
%% 2.汇总视频特征
Festure_time_freq = [  fest_five_band];
Festure_time_freq_candidate_num = size(Festure_time_freq,2);

end

