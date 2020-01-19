
function [realT, distr, prob] = permutationTwoGroups(vec1,vec2,numPerm)
% for debugging:
% vec1=same_gr1_all_subs_resh(1,:);
% vec2=same_gr2_all_subs_resh(1,:);
% numPerm=1000;
assert(length(vec1)==length(vec2),'vectos have to be the same length')
combVec=[vec1;vec2]';
[h, p, ci, stats]=ttest(vec1,vec2);
realT=stats.tstat;
distr=zeros(numPerm,1);
distr(1)=realT;
for iPerm=2:numPerm
    tempRandVec=zeros(size(combVec));
    for iVal=1:length(vec1)
        tempRand=randperm(2);
        tempRandVec(iVal,:)=combVec(iVal,tempRand);
    end
    [h, p, ci, stats]=ttest(tempRandVec(:,1),tempRandVec(:,2));
    distr(iPerm)=stats.tstat;
end
prob=sum(distr>realT)/numPerm;
end