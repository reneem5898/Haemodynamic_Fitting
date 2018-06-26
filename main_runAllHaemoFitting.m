%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script runs all haemo fitting for all pig datasets in one go
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Date modified: 08/02/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
diary off

parentDir = 'P:\Data\OSU Models\HFPEF';

% Number of frames in cine dataset
numFrames = 30;

% Use log file to save details of haemodynamic fitting
logFile = sprintf('%s/haemo_log_file_all.txt', parentDir);
delete(logFile); % Delete to start log file from scratch
diary(logFile);


% File with cardiac events chosen from MRI cine dataset
fid = fopen(sprintf('%s/OSU_Pig_Names_Cardiac_Events_all.txt', parentDir));
allLines = textscan(fid, '%s %s %d %d %d %d %d');


% Create matrix of relevant study frame info for material parameter optimisation
studyNamesFile = fopen(sprintf('%s/StudyNames.txt', parentDir), 'w');

% Initialise variables for pressures and volumes
p = zeros(size(allLines{1},1), 2, numFrames); % (# animals x 3 time points (bx, m1, m2)), # cardiac time points
v = zeros(size(allLines{1},1), numFrames); % (# animals, 3 time points (bx, m1, m2), # cardiac time points

% Initialise variables to save time points and number of frames
allTimePoints = {};
allStudyFrames = zeros(size(allLines{1},1),1);

for i = 1:size(allLines{1},1)
    
    close all
    
    fprintf('\n\n------------------------------------------------------------\n')
       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Cardiac events from CIM
    edMRI = allLines{3}(i);
    eivcMRI = allLines{4}(i);
    esMRI = allLines{5}(i);
    eivrMRI = allLines{6}(i);
    dsMRI = allLines{7}(i);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Construct folder names - specific to Ohio data and my naming convention
    % - CHANGE AS NEEDED
    
    % Animal name and time point (e.g. Pig1_ET38 and Base)
    animal = char(allLines{1}(i));
    tp = char(allLines{2}(i));
    disp(animal)
    disp(tp)
    
    % Directory to find pressure data and also where figures, mat files, etc. will be saved
    directory = sprintf('P:/Data/OSU Models/HFPEF/%s/%s/Young', animal, tp);
    if ~exist(directory, 'dir')
        directory = sprintf('P:/Data/OSU Models/HFPEF/%s/%s/Pressure Measurements', animal, tp);
    end
    
    % CIM model directory - where to find CIM model files for specific case/volunteer/animal/etc.
    model_directory = sprintf('C:/AMRG/CIM_81_WARP/CIM_DATA/%s_HFPEF_%s', animal, tp);
    if ~exist(model_directory, 'dir')
        model_directory = sprintf('C:/AMRG/CIM_81_WARP/CIM_DATA/%s_HFPEF_%s+mre', animal, tp);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Construct output pressure filename - for saving registered pressure
    % - FOR EASE OF USE WITH C1 PARAMETER ESTIMATION - Pressure data file should be: <model name>_registered_LVP.txt
    
    tmp = strsplit(model_directory, '/');
    outPressureFile = sprintf('%s/%s_registered_LVP.txt', directory, tmp{end});
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run haemodynamic fitting code - main function
    [pressure, volume, pressure_unshifted] = main_HaemoFitting(directory, model_directory, dsMRI, edMRI, esMRI, eivcMRI, eivrMRI, outPressureFile);
    p(i,1,1:size(pressure,1)) = pressure(:,2);
    p(i,2,1:size(pressure,1)) = pressure_unshifted;
    v(i,1:size(volume,1)) = volume;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save list of study names: <animal>_<tp>
    studyName = sprintf('%s', tmp{end});
    
    % Save DS, ED, ES and total number of frames
    studyFrames = [dsMRI, edMRI, esMRI, size(pressure,1)];
    
    % Create a full list of all time points (Bx, M1 and M2) and frame 
    % numbers for each study - used for organising plots at the end
    allTimePoints{i} = tp;
    allStudyFrames(i) = size(pressure,1);
    
    % Write study name and relevant frame numbers to a file
    fprintf(studyNamesFile, '%s\t', studyName);
    fprintf(studyNamesFile, '%d\t%d\t%d\t%d\n', studyFrames);
    
end

% Close study frames file
fclose(studyNamesFile);

% Create plots for shifted and unshifted pressures
plot_PressureVolumes(squeeze(p(:,1,:)), v, parentDir, [], [], allTimePoints, allStudyFrames);
plot_PressureVolumes(squeeze(p(:,2,:)), v, parentDir, '-Unshifted', [36], allTimePoints, allStudyFrames);

% Create plots - distinguishing results from each study
plot_PressureVolumes_Study(squeeze(p(:,1,:)), v, parentDir, [], [], allTimePoints, allStudyFrames);

% Create plots for shifted and unshifted pressures with cases removed which
% had low identifiability
% plot_PressureVolumes(squeeze(p(:,1,:)), v, parentDir, [], [1, 2, 3, 21], allTimePoints, allStudyFrames);
% plot_PressureVolumes(squeeze(p(:,2,:)), v, parentDir, '-Unshifted', [1, 2, 3, 21], allTimePoints, allStudyFrames);

% Save all resulting pressures and volumes
save(sprintf('%s/pressures-volumes.mat', parentDir), 'p', 'v', 'allTimePoints', 'allStudyFrames');

diary off