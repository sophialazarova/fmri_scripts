% path
path = 'G:\Analyse\Stress\Data\Experiment_fMRI_Jun2017-Okt2018'; 

% load data 
free_recall = xlsread([path '/FreeRecallMatlab.xlsx']); 
% 1: TSST or fTSST
% 2: central
% 3: peripheral
% 4: exclusion

% how many permutations and subjects?
permut = 1000; 
subject_nr = sum(free_recall(:,4) == 0); 

free_recall_excl = free_recall(free_recall(:,4) == 0,2:3); 
BS_excl = free_recall(free_recall(:,4) == 0,1); 

% actual anova 
[table,~] = simple_mixed_anova(free_recall_excl,BS_excl); 

% write down actual pvalues
Fval_BS = table2array(table(2,4)); 
Fval_WS = table2array(table(4,4)); 
Fval_int = table2array(table(5,4)); 

% prealloc
Fval_BS_nulldis = NaN(1,1000); 
Fval_WS_nulldis = NaN(1,1000); 
Fval_int_nulldis = NaN(1,1000); 

% distribution of p-values for shuffled data
for per = 1:permut
    
    disp(num2str(per)); 
    
    % create shuffled values 
    free_recall_to_shuffle = free_recall_excl; 
    free_recall_long = free_recall_to_shuffle(:); 
    free_recall_long_shuffled = free_recall_long(randperm(length(free_recall_long))); 
    free_recall_shuffled = reshape(free_recall_long_shuffled,subject_nr,2); 
    
    % carry out anova on shuffled values
    [table,~] = simple_mixed_anova(free_recall_shuffled,BS_excl); 
    
    % write down shuffled pvalues
    Fval_BS_nulldis(per) = table2array(table(2,4)); 
    Fval_WS_nulldis(per) = table2array(table(4,4)); 
    Fval_int_nulldis(per) = table2array(table(5,4)); 
    
end 

quantile(Fval_BS_nulldis,0.95)
quantile(Fval_WS_nulldis,0.95)
quantile(Fval_int_nulldis,0.95)

figure()
subplot(1,3,1)
histogram(Fval_WS_nulldis)
hold on; 
line([Fval_WS Fval_WS], [0 600],'Color','red');
title('Central vs. peripheral')
subplot(1,3,2)
histogram(Fval_BS_nulldis)
hold on; 
line([Fval_BS Fval_BS], [0 600],'Color','red');
title('Stress vs. Friendly')
subplot(1,3,3)
histogram(Fval_int_nulldis)
hold on; 
line([Fval_int Fval_int], [0 600],'Color','red');
title('Interaction')

