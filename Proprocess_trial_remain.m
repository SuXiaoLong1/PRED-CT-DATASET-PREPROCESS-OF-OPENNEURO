function remain_trial=Proprocess_trial_remain(Proprocess_trial_abandon_target)
%获取每个通道剩余数据量
%分类
remain_trial = zeros(1,size(Proprocess_trial_abandon_target,2));
for i=1:size(Proprocess_trial_abandon_target,2)
   for j=1:size(Proprocess_trial_abandon_target,1)
   if isempty(Proprocess_trial_abandon_target{j,i})==0
        remain_trial(i)=remain_trial(i)+1;
        
   end
   end
end
end

