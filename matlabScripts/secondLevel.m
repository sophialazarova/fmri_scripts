
%% second-level contrasts -------------------------------------------------
% 
% % change this with more subjects!
% participant_with_missing_condition = zeros(1,length(subjects));
% participant_with_missing_condition(2) = 1; 

%% a) Contrast all subs ---------------------------------------------------
% Create batch for model estimation
%prior to this step, all contrast files ( spmT_0001) are renamed
%(spmT_0001_subjName) and copied in another folder altogether. 
load(fullfile(currDir,'Stats','11_contrasts_one_sided_ttest_second_level_cogSci')) % load batch
matlabbatch{1}.spm.stats.factorial_design.dir{1} = fullfile(dataDir,'Analysis','secondLevel','con1');
count = 1; 
for sub = 1:length(subjects)

        % display how far you are
        fprintf('Subject: %s\n', subjects{sub, 1});
        
        % define contrasts
         matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{count} = fullfile(dataDir,'Analysis','secondLevel','con1',['spmT_0001_' subjects{sub,1} '.nii']);
%         matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{count} = fullfile(dataDir,'Analysis','secondLevel','con2',['spmT_0002_' subjects{sub,1} '.nii']);
%          matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{count} = fullfile(dataDir,'Analysis','secondLevel','con3',['spmT_0003_' subjects{sub,1} '.nii']);
%         matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{count} = fullfile(dataDir,'Analysis','secondLevel','con4',['spmT_0004_' subjects{sub,1} '.nii']);

%         matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';  %replication of contrast for other sessions
        count=count+1;
end
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = transpose(matlabbatch{1}.spm.stats.factorial_design.des.t1.scans);

% execute batch
spm_jobman('initcfg')
spm_jobman('run',matlabbatch)
clear matlabbatch

%% second-level estimation ------------------------------------------------
%% at the moment this is done manually
% Create batch for model estimation
load(fullfile(currDir,'Stats','12_second_level_est')) % load batch
matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile(dataDir,'Analysis','secondLevel','con1','SPM.mat');
spm_jobman('initcfg')
spm_jobman('run',matlabbatch)
clear matlabbatch
load(fullfile(currDir,'Stats','12_second_level_est')) % load batch
matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile(dataDir,'Analysis','secondLevel','con2','SPM.mat');
spm_jobman('initcfg')
spm_jobman('run',matlabbatch)
clear matlabbatch
load(fullfile(currDir,'Stats','12_second_level_est')) % load batch
matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile(dataDir,'Analysis','secondLevel','con3','SPM.mat');
spm_jobman('initcfg')
spm_jobman('run',matlabbatch)
clear matlabbatch
load(fullfile(currDir,'Stats','12_second_level_est')) % load batch
matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile(dataDir,'Analysis','secondLevel','con4','SPM.mat');

% execute batch
spm_jobman('initcfg')
spm_jobman('run',matlabbatch)
clear matlabbatch


% %% copy files
% for sub = 1:length(subjects)
%     mkdir(fullfile('I:\Daten_fuer_Annika',filesep,subjects{sub,1}))
%     copyfile(fullfile(main_path,filesep,'fMRIData',filesep,subjects{sub,1},filesep,'GLM_uni'),fullfile('I:\Daten_fuer_Annika',filesep,subjects{sub,1}))
% end


% do by hand! 
%% ROI-Analysis -----------------------------------------------------------
% take ROIs from wfu-pickatlas and save mean activity of the ROIs 
ROI_names = {'Fusiform','Amygdala','PPA'}; 
for sub = 1:length(subjects)
    
    % display how far you are
    fprintf('Subject: %s\n', subjects{sub, 1});
    
    for iROI = 1 %1:length(ROI_name)
        
        % how would you like to call your output-file?
        ofilename = ROI_names{iROI};
        
        % get the data of the specific ROI file
        V_ROI = spm_vol(strcat(currDir, filesep, 'masks', filesep, ['r' ofilename '_ROI.nii']));
        [Y_ROI, ~] = spm_read_vols(V_ROI); % Y = data, XYZ = coordinates
        
        % which index corresponds to the hippocampus
        idx_ROI = [];
        idx_ROI = [idx_ROI ; find(Y_ROI == 1)];
        
        % get the data from the contrast file 
        V_t_sub = spm_vol(strcat(dataDir, filesep, 'Analysis', filesep, 'secondLevel', filesep, 'spmT_0001.nii'));
        [Y_t_sub, ~] = spm_read_vols(V_t_sub); % Y = data, XYZ = coordinates
        
%         % check how the images fit together
        plotImageAndOverlayColor(Y_t_sub,Y_ROI,'fig',46,3,[1 1 1],[46 46 46])
        
        % extract the values out of the mask from the t-map
        ROI_values = Y_t_sub(idx_ROI);
        meanvaluesperROI(iROI,sub) = mean(ROI_values);
        
    end
end