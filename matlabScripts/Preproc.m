%% fMRI MarsPrimefMRI -------------------------------------------------
% Subjects viewed pictures of the martian objects. Their
% task was to find the colored fixation cross inbetween, but they were instructed
% to look at the objects carefully. During each run,each object (12
% objects) was presented once in random order. There were 10 runs.
% Each run was interrupted by a 25 second break in which pictures of houses
% or pictures of faces were presented (task for fMRI practical)

% Get information about data quality: it might be good to see how many red
% crosses were missed - TASK FOR LATER (Info about sleep on questionnaires)

% Preprocessing:
% - Slice time correction (ascending): I checked and I am sure it was an
% ascending slice acquisition protocol (see script:
% slice_acquisition_check)
% - Motion correction: Subjects/Runs are excluded if subjects move more
% than one voxel-size (2.5mm) during a run and/or if the motion profile is
% very spiky). This information is again written in subjects-variable.
% - Smooting: Just in case you would like to try the MVPA-analysis with
% minimally smoothed data (kernel: 5mm)

% GLM:
% They fly and each object is modeled as a regressor (73 conditions per
% run). The first derivative of each regressor is also modeled (73
% conditions per run). In addition, we use the motion parameters as
% regressors and the mean of each run (7 regressors). Overall, there are 6
% more constant regressors.

% Create ROIs from aparc-aseg template:
% The specified ROIs are derived from the aparc-aseg template of each
% participant and coregistered to the meanEPI in order to extract the
% voxel-time-courses of those ROIs.

% MVPA analysis/RSA analysis:
% Extract Beta-values for each condition and ROI (transform to t-values by
% dividing by MS_res). Correlate each Beta-value with every other
% beta_value. Take only correlations between runs and calculate the mean
% value of each condition. Plot the RSA-matrices. Calculate averages for
% the conditions (central vs. peripheral objects) and for memory conditions
% (hits vs. misses) of on- and off-diagonal correlations values (same vs.
% diff).

%% paths and variables ----------------------------------------------------
clear
clc
cd /Users/sophia/Desktop/matlabScripts
% set path
currDir=pwd;
currDirSplit=strsplit(currDir,filesep);
if isunix
    baseDir=['/' fullfile(currDirSplit{1:end-1})];
elseif ispc
    baseDir=fullfile(currDirSplit{1:end-1});
%elseif ismac
%    baseDir=fullfile(currDirSplit{1:end-1});
%    warning('Filenaming might not work on MAC')
end
display(['Base directory is ' baseDir])

dataDir='/Volumes/ExternalLoo/MRIPracticalData';
behavDirSL='/Volumes/ExternalLoo/MRIPracticalData/beh_fmriTask/shortLearn';
behavDirLL='/Volumes/ExternalLoo/MRIPracticalData/beh_fmriTask/longLearn';
dicomDir=fullfile(dataDir,'DICOM'); % no need
nifti3DDir='/Volumes/ExternalLoo/MRIPracticalData/mriConvert';
nifti4DDir=fullfile(dataDir,'nifti'); % no need

SPMpath='/Users/sophia/software/spm12';
addpath(SPMpath);
stDir=fullfile(dataDir,'stcorr');

subjects = {
    '001'
    '002'
    '003'
    '004'
    '005'
    '006'
    '007'
    '008'
    '009'
    '010'
    '011'
    '012'
    '013'
    '014'
    '051'
    '052'
    '053'
    '054'
    '055'
    '056'
    '057'
    '058'
    '059'
    '060'
    '061'
    '062'
    '063'
    '064'
    '065'
    '066'
    '067'
    '068'
    '069'
    '070'
     };



% subjects
% es darf niemals ein Proband auskommentiert werden!! (sonst verschieben
% sich die Probanden nach unten und die behavioralen Daten stimmen nicht
% mehr)


% define ROIs
label_value_all = {
    [1000:2999] % whole cortex
    %     [17,53],... % hippocampus
    %     [17],... % hippocampus left
    %     [53],... & hippocampus right
    %     [18,54],... % amygdala
    %     [18],... % amygdala left
    %     [54],... % amygdala right
    [1021,2021] % pericalcarine
    [1007,2007] % fusiform
    %     [1013,2013],... % lingual
    %     [1011,2011],... % lateraloccipital
    %     [1005,2005],... % cuneus
    %     [1025,2025],... % precuneus
    %     [1030,2030],... % superiortemporal
    %     [1015,2015],... % middletemporal
    %     [1009,2009],... % inferiortemporal
    [1016,2016] % parahippocampal
    %     [1033,2033],... % temporalpole
    %     [1016,2016,1021,2021,1005,2005,1007,2007,1011,2011,1013,2013]}; % VVS:
    % parahippocampal, pericalcarine, cuneus, fusiform, lateraloccipital, lingual
    };

ROI_name_all = {
    %     'whole_GM',...
    'whole_cortex'
    %     'hippocampus',...
    % 'hippocampus_anterior',...
    % 'hippocampus_middle',...
    % 'hippocampus_posterior',...
    % 'hippocampus_anterior_left',...
    % 'hippocampus_anterior_right',...
    % 'hippocampus_middle_left',...
    % 'hippocampus_middle_right',...
    % 'hippocampus_posterior_left',...
    % 'hippocampus_posterior_right',...
    %     'hippocampus_left',...
    %     'hippocampus_right',...
    %     'amygdala',...
    %     'amygdala_left',...
    %     'amygdala_right',...
    'pericalcarine'
    'fusiform'
    %     'lingual',...
    %     'lateraloccipital',...
    %     'cuneus',...
    %     'precuneus',...
    %     'superiortemporal',...
    %     'middletemporal',...
    %     'inferiortemporal',...
    'parahippocampal'
    %     'temporalpole',...
    %     'VVS',...
    };
%
colors = ['r', 'b', 'g', 'c', 'm', 'k'];

% % todays_date = datetime(t2, 'ConvertForm','datenum');

%% preprocessing ----------------------------------------------------------
%
%% slice time correction
for sub = 1:size(subjects,1)

    % display how far you are
    fprintf('Subject: %s\n', subjects{sub, 1});
    participantFolder=fullfile(nifti3DDir,subjects{sub,1});

    subFolders=dir(participantFolder);

    % specify fMRI data
     fils = dir(fullfile(participantFolder,filesep,['MarsPrimefMRI_' subjects{sub} '_1'])); fils = fils(~[fils.isdir]);  % save data path into "fils"


    %% a. slice time correction (batch)
    load(strcat(currDir,filesep,'Preprocessing',filesep,'3_slice_time_correction_ascending')) % load batch

     matlabbatch{1}.spm.temporal.st.scans(:) = [];

    % put data into batch
    % scans
    count=1;
    for scan = 1:length(fils)
        matlabbatch{1}.spm.temporal.st.scans{count}{scan,1} = fullfile(participantFolder,['MarsPrimefMRI_' subjects{sub} '_1'],fils(scan).name);  %participantFolder,['MarsPrimefMRI_' subjects{sub} '_1'],
    end

    % execute batch
    spm_jobman('initcfg')
    spm_jobman('run',matlabbatch)

    % create folders and move data into folders
    mkdir(strcat(dataDir,filesep,'stCorr',filesep,subjects{sub,1}));
    movefile(strcat(participantFolder,filesep,['MarsPrimefMRI_' subjects{sub,1} '_1'],filesep,'a*'),strcat(dataDir,filesep,'stcorr',filesep,subjects{sub}));
    clear matlabbatch fils;
end


%% motion correction
for sub = 1:size(subjects,1)
    
    % display how far you are
    fprintf('Subject: %s\n', subjects{sub, 1});
    
    % specify fMRI data
    fils = dir(fullfile(dataDir,filesep,'stCorr',filesep,subjects{sub})); fils = fils(~[fils.isdir]); % save data path into "fils"
    
    %% r. motion correction (batch)
    load(strcat(currDir,filesep,'Preprocessing',filesep,'4_motion_correction')) % load batch
    
    matlabbatch{1}.spm.spatial.realign.estwrite.data(:) = [];
    % put data into batch
    % scans
    count = 1;
    for scan = 1:length(fils)
        matlabbatch{1}.spm.spatial.realign.estwrite.data{count}{scan,1} = fullfile(dataDir,'stcorr',subjects{sub},fils(scan).name);
    end
    
    % execute batch
    spm_jobman('initcfg')
    spm_jobman('run',matlabbatch)
    
    % create folders and move data into folders
    mkdir(strcat(dataDir,filesep,'moCorr',filesep,subjects{sub,1}));
    movefile(strcat(dataDir,filesep,'stCorr',filesep,subjects{sub},filesep,'ra*.nii'),strcat(dataDir,filesep,'moCorr',filesep,subjects{sub,1}));
    movefile(strcat(dataDir,filesep,'stCorr',filesep,subjects{sub},filesep,'rp*.txt'),strcat(dataDir,filesep,'moCorr',filesep,subjects{sub,1}));
    movefile(strcat(dataDir,filesep,'stCorr',filesep,subjects{sub},filesep,'mean*.nii'),strcat(dataDir,filesep,'moCorr',filesep,subjects{sub,1}));
    
    % take a look at the motion parameters: load parameter files
    motion_fils = dir(fullfile(dataDir,'moCorr',subjects{sub},'rp_*')); motion_fils = motion_fils(~[motion_fils.isdir]); % save data path into "motion_fils"
    motion_params = load(fullfile(dataDir,'moCorr',subjects{sub},motion_fils.name));
    
    
    % create a figure
    h = figure;
    h.Name = subjects{sub};
    count = 1;
    
    for param = 1:6
        plot(motion_params(:,param),colors(param));
        hold on;
    end
    ylim([-3 3])
    legend({'trans x', 'trans y', 'trans z', 'pitch', 'roll', 'yaw'},'Location','southeastoutside');
    title([subjects{sub}])
    
    
    saveas(h, fullfile(dataDir,'moCorr',subjects{sub},[subjects{sub} '_motion_params.bmp']))
    close;
    
    clear matlabbatch fils;
end


% sub = 1
%% Coregistration
for sub = 1:size(subjects,1)
    
    % display how far you are
    fprintf('Subject: %s\n', subjects{sub, 1});
    
    %Load scans:
    %cd(strcat(mainpath,num2str(i),'\ra'));
    fils = dir(fullfile(dataDir,filesep,'moCorr',filesep,subjects{sub},filesep,'ra*.nii')); fils = fils(~[fils.isdir]); % save data path into "fils"
    % ra_scans=dir('*.nii');
    ra_names={fils.name};
    
    for j=1:length(ra_names)
        ra_files(j,1)=strcat(dataDir,'\moCorr\',subjects{sub},filesep,ra_names(j));
    end
    
    path_meanEPI = fullfile(dataDir,filesep,'moCorr',filesep,subjects{sub});
    
    % cd(strcat(mainpath,num2str(i),'\T1'));
    T1_scan=dir(fullfile(nifti3DDir,filesep,subjects{sub},filesep,'*_Anat_*',filesep,'*.nii'));
    T1_name={T1_scan.name};
    T1=fullfile(T1_scan.folder,T1_name);
    %strcat(nifti3DDir,filesep,subjects{sub},filesep,'*_Anat_*',filesep,'*_Anat_*.nii');
    
    % cd(strcat(dataDir,num2str(i),'\a'));
    Mean_scan=dir(fullfile(dataDir,filesep,'moCorr',filesep,subjects{sub},filesep,'mean*.nii'));
    Mean_name={Mean_scan.name};
    Mean_file=strcat(dataDir,'\moCorr',filesep,subjects{sub},filesep,Mean_name);
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = T1;
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = Mean_file;
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = ra_files;
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 7;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    
    %Load and start the batch:
    %save ('matlabbatch');
    % spm_jobman('initcfg')
    spm_jobman('run', matlabbatch)
    
    mkdir(strcat(dataDir,filesep,'coreg',filesep,subjects{sub,1}));
    movefile(strcat(dataDir,filesep,'moCorr',filesep,subjects{sub},filesep,'rra*.nii'),strcat(dataDir,filesep,'coreg',filesep,subjects{sub,1}));
    
    clear matlabbatch fils;
    
end

%% Normalization
for sub = 1:size(subjects,1)
    % display how far you are
    fprintf('Subject: %s\n', subjects{sub, 1});
    
    
    %Load scans:
    
    fils = dir(fullfile(dataDir,filesep,'coreg',filesep,subjects{sub},filesep,'rra*.nii')); fils = fils(~[fils.isdir]); % save data path into "fils"
    rra_names={fils.name};
    
    for j=1:length(rra_names)
        rra_files(j,1)=strcat(dataDir,'\coreg\',subjects{sub},filesep,rra_names(j));
    end
    
    T1_scan=dir(fullfile(nifti3DDir,filesep,subjects{sub},filesep,'*_Anat_*',filesep,'*.nii'));
    T1_name={T1_scan.name};
    T1=fullfile(T1_scan.folder,T1_name);
    
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = T1;
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = rra_files;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {strcat(SPMpath,'spm12\tpm\TPM.nii')};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 7;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
    
    %Load and start the batch:
    %save ('matlabbatch');
    spm_jobman('run', matlabbatch)
    
    mkdir(strcat(dataDir,filesep,'norm',filesep,subjects{sub,1}));
    movefile(strcat(dataDir,filesep,'coreg',filesep,subjects{sub},filesep,'wrra*.nii'),strcat(dataDir,filesep,'norm',filesep,subjects{sub,1}));
    
    clear matlabbatch fils;
    
end

%% Smoothing
for sub = 1:size(subjects,1)
    %Load scans:
    fils = dir(fullfile(dataDir,filesep,'norm',filesep,subjects{sub},filesep,'wrra*.nii')); fils = fils(~[fils.isdir]); % save data path into "fils"
    wrra_names={fils.name};
    
    for j=1:length(wrra_names)
        wrra_files(j,1)=strcat(dataDir,'\norm\',subjects{sub},filesep,wrra_names(j));
    end
    
    %wrra_files=[wrraFilm_files; wrraRest_files; wrraRT_files];
    
    matlabbatch{1}.spm.spatial.smooth.data = wrra_files;
    matlabbatch{1}.spm.spatial.smooth.fwhm = [5 5 5];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    
    %Load and start the batch:
    % save ('matlabbatch');
    spm_jobman('run', matlabbatch)
    
    mkdir(strcat(dataDir,filesep,'smo',filesep,subjects{sub,1}));
    movefile(strcat(dataDir,filesep,'norm',filesep,subjects{sub},filesep,'swrra*.nii'),strcat(dataDir,filesep,'smo',filesep,subjects{sub,1}));
    clear matlabbatch fils;
    
end
