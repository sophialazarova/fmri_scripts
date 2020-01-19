function [realT, distr, prob] = permutationObjects(same_per_object,diff_per_object,objects,numPerm)
% for debugging:
% same_per_object=same_per_object_all_subs(:,ind_mask,:);
% diff_per_object=diff_per_object_all_subs(:,:,ind_mask,:);
% numPerm=1000;
% objects = objects; 

% real T 
group1 = objects(1:9);
group2 = objects(10:end); 
same_diff_group1_mean = mean(same_per_object(group1,:),1)' - squeeze(mean(nanmean(diff_per_object(group1,:,:),2),1));
same_diff_group2_mean = mean(same_per_object(group2,:),1)' - squeeze(mean(nanmean(diff_per_object(group2,:,:),2),1));
[h, p, ci, stats]=ttest(same_diff_group1_mean,same_diff_group2_mean);
realT=stats.tstat;

% distribution of Ts
distr=zeros(numPerm,1);
distr(1)=realT;
for iPerm=2:numPerm
    objects = objects(randperm(length(objects)));
    group1 = objects(1:9);
    group2 = objects(10:end); 
    same_diff_group1_mean = mean(same_per_object(group1,:),1)' - squeeze(mean(nanmean(diff_per_object(group1,:,:),2),1));
    same_diff_group2_mean = mean(same_per_object(group2,:),1)' - squeeze(mean(nanmean(diff_per_object(group2,:,:),2),1));
    [h, p, ci, stats]=ttest(same_diff_group1_mean,same_diff_group2_mean);
    distr(iPerm)=stats.tstat;
end
prob=sum(distr>realT)/numPerm;
end