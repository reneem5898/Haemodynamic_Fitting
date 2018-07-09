function plot_PressureVolumes_Study(p, v, parentDir, descriptor, StudyFrames)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots a) pressure-volume loops, b) plots from diastasis (DS)
% to end-diastole (ED) and c) plots of only DS and ED
% 
% Inputs: 1) p - pressure (kPa)
%         2) v - volume (mL)
%         3) parentDir - directory where to save plots
%         4) descriptor - string to attach to end of files (e.g. "Study1",
%         "Unshifted", etc.)
%         7) StudyFrames - vector of integers indicating the number of
%         frames in each study
%
% Plots created: 1) pressure-volume loops
%                2) diastolic pressure curve (from DS to ED)
%                3) pressure at DS and ED only
%
% Renee Miller (renee.miller@auckland.ac.nz)
% Date Modified: 25 June 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of cases
numCases = size(p, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot all pressure-volume loops together
FH = figure;
for i = 1:numCases
    
    % Get pressure and volume for current study
    volume = v(i,1:StudyFrames(i,3));
    pressure = p(i,1:StudyFrames(i,3));
    
    if StudyFrames(i,3) == 30
        h1 = plot(volume, pressure, 'm-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
    else
        h2 = plot(volume, pressure, 'k-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end
    hold on
end
xlabel('Volume (mL)', 'FontSize', 12)
ylabel('Pressure (kPa)', 'FontSize', 12)
legend([h1, h2], {'Study 1', 'Study 2'}, 'Location', 'Best')
saveas(FH, sprintf('%s/PV-Loops-Study%s.png', parentDir, descriptor));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot all diastolic portions of the PV loops together
FH2 = figure;
for i = 1:numCases
    
    % Get pressure and volume for current study
    volume = v(i,1:StudyFrames(i,3));
    pressure = p(i,1:StudyFrames(i,3));
    
    % Get diastasis and end diastolic time points
    dsMRI = StudyFrames(i,1);
    edMRI = StudyFrames(i,2);
    
    % Get only volume and pressure during diastole: from DS to ED
    % (inclusive)
    if dsMRI > edMRI
        volume_diastole = [volume(dsMRI:end), volume(1:edMRI)];
        pressure_diastole = [pressure(dsMRI:end), pressure(1:edMRI)];
    else
        volume_diastole = volume(dsMRI:edMRI);
        pressure_diastole = pressure(dsMRI:edMRI);
    end
    
    if StudyFrames(i,3) == 30
        h1 = plot(volume_diastole, pressure_diastole, 'm-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
    else
        h2 = plot(volume_diastole, pressure_diastole, 'k-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end
    hold on
    
end
xlabel('Volume (mL)', 'FontSize', 12)
ylabel('Pressure (kPa)', 'FontSize', 12)
legend([h1, h2], {'Study 1', 'Study 2'}, 'Location', 'Best')
saveas(FH2, sprintf('%s/PV-Diastole-Study%s.png', parentDir, descriptor));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot DS and ED
FH3 = figure;
for i = 1:numCases
    
    % Get pressure and volume for current study
    volume = v(i,1:StudyFrames(i,3));
    pressure = p(i,1:StudyFrames(i,3));
    
    % Get diastasis and end diastolic time points
    dsMRI = StudyFrames(i,1);
    edMRI = StudyFrames(i,2);
    
    % Get only volume and pressure during diastole: from DS to ED
    % (inclusive)
    if dsMRI > edMRI
        volume_diastole = [volume(dsMRI:end), volume(1:edMRI)];
        pressure_diastole = [pressure(dsMRI:end), pressure(1:edMRI)];
    else
        volume_diastole = volume(dsMRI:edMRI);
        pressure_diastole = pressure(dsMRI:edMRI);
    end
    
    if StudyFrames(i,3) == 30
        h1 = plot(volume_diastole([1,end]), pressure_diastole([1,end]), 'm-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
    else
        h2 = plot(volume_diastole([1,end]), pressure_diastole([1,end]), 'k-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end
    hold on
    
end

xlabel('Volume (mL)', 'FontSize', 12)
ylabel('Pressure (kPa)', 'FontSize', 12)
legend([h1, h2], {'Study 1', 'Study 2'}, 'Location', 'Best')
saveas(FH3, sprintf('%s/PV-Diastole-DS-ED-Study%s.png', parentDir, descriptor));

close all