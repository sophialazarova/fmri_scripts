
%% GLM (specification)-----------------------------------------------------
for sub = 1:size(subjects,1)

    fils = ls(fullfile(behavDirSL,filesep,'*_',subjects{sub},'_*')); % save data path into "fils"
    fils = [];

    % display how far you are
    disp(['Subject: ', subjects{sub, 1}]);

    % read behavioral data - run once for long learn and then for short
    % learn!!
     beh_path = strcat(behavDirSL,filesep,'fmripractical_shortLearn_',subjects{sub}, '_1.txt'); %beh_path = beh_path(~[beh_path.isdir]);
    %beh_path = strcat(behavDirLL,filesep,'fmripractical_longLearn_vp_',subjects{sub}, '_1.txt'); %beh_path = beh_path(~[beh_path.isdir]);
    fileID = fopen(beh_path);
    beh_data(:,:) = textscan(fileID, '%d %d %d %d %d %d %d %d %d %d %d %d %d %d', 'HeaderLines',1);  % header: trial_kind	trial_no	pulse_no	firstFixStart	trial_start_time_pulse_abs	picture_presented	category_of_picture	secondFixStart	fix_colored	redfix_onset_time	response	response_time	1or2TR	thirdFixTimeStart	
    
    % read fMRI data
    fils = dir(fullfile(dataDir,'smo',subjects{sub,1}, 'swrra*')); fils = fils(~[fils.isdir]); % save data path into "fils"
    
    % create batch for model specification
    load(strcat(currDir,filesep,'GLM',filesep,'8_first_level_specs')) % load batch
    
    % fill in the batch
    
    % fill in directory
    matlabbatch{1}.spm.stats.fmri_spec.dir{1} = strcat(dataDir,filesep,'smo'); % load directory into batch
    
    % fill in sess
    count = 1;

    % scans
    for scan = 1:length(fils)
        matlabbatch{1}.spm.stats.fmri_spec.sess(count).scans{scan} = fullfile(dataDir,'smo',subjects{sub},fils(scan).name);
    end

    % conditions
    % every pic is one condition
    conditions = [1; 2; 3]; % 
    numOnset= 20;

      beh_data{1,15} = beh_data{1,4}+5060;
    
    
    trialKindMat=cell2mat(beh_data(1,1));
    trialNumbMat=cell2mat(beh_data(1,2));
    stimTimeMat=cell2mat(beh_data(1,15));
    
    
    onsetTime = zeros(3,10);
    durationTime = zeros(3,10);
    onsetsMars=(1:12:120);
    onsetMarsCount=1;
    onsetsBreak = (1:10:50);
    onsetFaceCount=1;
    onsetBuildCount=1;
        
    %defining 3 separate blocks - one for faces, one for buildings and one for
    %martian objects (detection ?)
    for trial = 1:length(trialKindMat(:,1))
        if trialKindMat(trial,1)==1 && trialNumbMat(trial,1)==onsetsMars(1,onsetMarsCount)
                onsetTime(1,onsetMarsCount)=stimTimeMat(trial,1);
                if onsetMarsCount ==1
                   durationTime(1,onsetMarsCount)= stimTimeMat(trial+12,1) - stimTimeMat(trial,1) - 2500; 
                    onsetMarsCount=onsetMarsCount+1;                   
                elseif onsetMarsCount > 1 && onsetMarsCount < 10
                   durationTime(1,onsetMarsCount)= stimTimeMat(trial+12,1) - stimTimeMat(trial,1) - 2500;                                        
                    onsetMarsCount=onsetMarsCount+1;
                elseif onsetMarsCount == 10
                   durationTime(1,onsetMarsCount)= stimTimeMat(trial+12,1) - stimTimeMat(trial,1) - 2500;                                        
                    onsetMarsCount=1;
                end
        elseif trialKindMat(trial,1)==2 && trialNumbMat(trial,1)==onsetsBreak(1,onsetFaceCount)
                onsetTime(2,onsetFaceCount)=stimTimeMat(trial,1);
                
                if trial+10< length(trialKindMat(:,1))                
                    durationTime(2,onsetFaceCount)= stimTimeMat(trial+10,1) - stimTimeMat(trial,1) -2500;
                elseif trial+10 > length(trialKindMat(:,1))
                    durationTime(2,onsetFaceCount)= stimTimeMat(trial+9,1) + 2030 - stimTimeMat(trial,1);
                end 
                
                if onsetFaceCount < 5
                    onsetFaceCount=onsetFaceCount+1;
                elseif onsetFaceCount == 5
                    onsetFaceCount=1;
                end
        elseif trialKindMat(trial,1)==3 && trialNumbMat(trial,1)==onsetsBreak(1,onsetBuildCount)
                onsetTime(3,onsetBuildCount)=stimTimeMat(trial,1);
                
                if trial+10< length(trialKindMat(:,1))                
                    durationTime(3,onsetBuildCount)= stimTimeMat(trial+10,1) - stimTimeMat(trial,1) -2500;
                elseif trial+10 > length(trialKindMat(:,1))
                    durationTime(3,onsetBuildCount)= stimTimeMat(trial+9,1) + 2030 - stimTimeMat(trial,1);
                end
                
                if onsetBuildCount < 5
                    onsetBuildCount=onsetBuildCount+1;
                elseif onsetBuildCount == 5
                    onsetBuildCount=1;
                end
        end
    end
    
     % multiple condition file
    names = [];
    onsets = [];
    durations =[];
    
    for cond = 1:length(conditions)
        names{cond} = conditions(cond);
        onsets{cond} = onsetTime(cond,onsetTime(cond,:) > 0)/1000;
        durations{cond} = durationTime(cond,durationTime(cond,:) > 0)/1000;
    end

    % save multiple condition file
    save(fullfile(dataDir,filesep,'smo',filesep,subjects{sub},filesep,['MultiCond_' num2str(subjects{sub})]),...
        'names', 'onsets', 'durations')

    % fill in condition file
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi{1} = fullfile(dataDir,filesep,'smo',filesep,subjects{sub},filesep,['MultiCond_' num2str(subjects{sub}) '.mat']);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});

    % motion files (in order to specify nuisance regressors)
    motion_file = dir(fullfile(dataDir,filesep,'moCorr',filesep,subjects{sub,1},filesep,filesep,'rp*.txt')); 
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg{1} = fullfile(dataDir,filesep,'moCorr',filesep,subjects{sub},filesep,[motion_file.name]); % load motion parameters into batch           

    % execute batch
    spm_jobman('initcfg')
    spm_jobman('run',matlabbatch)

    % move spm.mat-file
    mkdir(fullfile(dataDir,filesep,'Analysis',filesep,'firstLevel',filesep,subjects{sub}));
%     delete(fullfile(main_path,filesep,'fMRIData',filesep,subjects{sub,1},filesep,'GLM',filesep,'*'));
    movefile(fullfile(dataDir,'smo','SPM.mat'),fullfile(dataDir,'Analysis', 'firstLevel',subjects{sub}))
    clear matlabbatch
    clear beh_data
end

%% GLM (estimation)--------------------------------------------------------
for sub = 1:size(subjects,1)

    % display how far you are
    fprintf('Subject: %s\n', subjects{sub, 1});

    % Create batch for model estimation
    load(fullfile(currDir,'GLM','9_first_level_est')) % load batch
    matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile(dataDir,'Analysis','firstLevel',subjects{sub},'SPM.mat');

    % execute batch
    spm_jobman('initcfg')
    spm_jobman('run',matlabbatch)

    % move beta-files etc. to a new folder
%     movefile(fullfile(dataDir,'smo',subjects{sub},'beta_*.nii'),fullfile(dataDir,'GLM',subjects{sub}))
%     movefile(fullfile(dataDir,filesep,'smo',filesep,subjects{sub},filesep,'ResMS.nii'),fullfile(dataDir,filesep,'GLM',filesep,subjects{sub}))
%     movefile(fullfile(dataDir,filesep,'smo',filesep,subjects{sub},filesep,'mask.nii'),fullfile(dataDir,filesep,'GLM',filesep,subjects{sub}))
%     movefile(fullfile(dataDir,filesep,'smo',filesep,subjects{sub},filesep,'RPV.nii'),fullfile(dataDir,filesep,'GLM',filesep,subjects{sub}))
%     movefile(fullfile(dataDir,filesep,'smo',filesep,subjects{sub},filesep,'SPM.mat'),fullfile(dataDir,filesep,'GLM',filesep,subjects{sub}))
%     
    clear matlabbatch

end


%% single-subject contrasts -----------------------------------------------
for sub = 1:length(subjects)
    
    % delete existing contrast files 
     delete(fullfile(dataDir,'Analysis','firstLevel',subjects{sub,1},'spmT_*.nii'));
     delete(fullfile(dataDir,'Analysis','firstLevel',subjects{sub,1},'con_*.nii'));
    
    % display how far you are
    fprintf('Subject: %s\n', subjects{sub, 1});
    
    % Create batch for model estimation
    load(fullfile(currDir,'Stats','10_contrasts_first_level')) % load batch
    matlabbatch{1}.spm.stats.con.spmmat{1} = fullfile(dataDir,'Analysis','firstLevel',subjects{sub,1},'SPM.mat');
    
    % define contrasts, maybe create another contrast, type 0 for all 
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Face>Build';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = repmat([0 1 -1 0 0 0],1);
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Build>Face';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = repmat([0 -1 1 0 0 0],1);
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Break>Mars';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = repmat([-1 0.5 0.5 0 0 0],1);
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Mars>Break';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = repmat([1 -0.5 -0.5 0 0 0],1);
                  
    % execute batch
    spm_jobman('initcfg')
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
end 
