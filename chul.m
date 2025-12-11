clear
clc
%%ftd----mat
for i =1:122
    file = ['sub-',sprintf('%03d', i),'_task-ProbabilisticSelection_eeg.set'];
    path = ['F:\python\EEG\data\ds003474-download\sub-',sprintf('%03d', i),'\eeg'];
    eeg = pop_loadset(file,path);
    data=eeg.data;
    save(['sub_',sprintf('%03d', i),'.mat'],'data');
end
clear
clc

