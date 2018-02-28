%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script runs all haemo fitting for all pig datasets in one go
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Date modified: 08/02/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
diary off

% Use log file to save details of haemodynamic fitting
logFile = 'P:/UA - PhD/OSU_Optimisation/haemo_log_file.txt';
delete(logFile); % Delete to start log file from scratch
diary(logFile);


% File with cardiac events chosen from MRI cine dataset
fid = fopen('P:/UA - PhD/OSU_Optimisation/OSU_Pig_Names_Cardiac_Events.txt');
allLines = textscan(fid, '%s %s %d %d %d %d %d');


% Create matrix of relevant study frame info for material parameter optimisation
studyNamesFile = fopen('P:/UA - PhD/OSU_Optimisation/StudyNames.txt', 'w');


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
    
    % Get subdirectories in image dir
    AllFiles = dir(directory);
    filenames = {AllFiles.name};
    subdirs = filenames([AllFiles.isdir]);
    subdirs(ismember(subdirs,{'.','..'})) = [];
    
    for j = 1:length(subdirs)
        if strfind(subdirs{j}, 'Cine_SAX')
            % Image directory - where to find MR SHORT AXIS images for dataset
            MR_directory = sprintf('%s/%s', directory, subdirs{j});
        end
    end
    
    % CIM model directory - where to find CIM model files for specific case/volunteer/animal/etc.
    model_directory = sprintf('C:/AMRG/CIM_81_WARP/CIM_DATA/%s_HFPEF_%s', animal, tp);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Construct output pressure filename - for saving registered pressure
    % - FOR EASE OF USE WITH C1 PARAMETER ESTIMATION - Pressure data file should be: <model name>_registered_LVP.txt
    
    tmp = strsplit(model_directory, '/');
    outPressureFile = sprintf('%s/%s_registered_LVP.txt', directory, tmp{end});
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run haemodynamic fitting code - main function
    [dsMRI, edMRI, esMRI] = main_HaemoFitting(directory, model_directory, dsMRI, edMRI, esMRI, eivcMRI, eivrMRI, outPressureFile);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save list of study names: <animal>_<tp>
    studyName = sprintf('%s', tmp{end});
    
    % Save DS, ED, ES and total number of frames
    studyFrames = [dsMRI, edMRI, esMRI, 30];
    
    % Write study name and relevant frame numbers to a file
    fprintf(studyNamesFile, '%s\t', studyName);
    fprintf(studyNamesFile, '%d\t%d\t%d\t%d\n', studyFrames);
    
end

% Close study frames file
fclose(studyNamesFile);

diary off