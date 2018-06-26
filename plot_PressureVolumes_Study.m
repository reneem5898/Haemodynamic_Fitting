function plot_PressureVolumes_Study(p, v, parentDir, descriptor, omit, TimePoints, StudyFrames)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots a) pressure-volume loops, b) plots from diastasis (DS)
% to end-diastole (ED) and c) plots of only DS and ED
% 
% Inputs: 1) p - pressure (kPa)
%         2) v - volume (mL)
%         3) parentDir - directory where to save plots
%         4) descriptor - string to attach to end of files (e.g. "Study1",
%         "Unshifted", etc.)
%         5) omit - studies to omit from plot if any (vector of numerical
%         values)
%         6) TimePoints - cell array of  values containing: 'Base', 'M1'
%         and 'M2' corresponding to the time point of each study
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
    volume = v(i,1:StudyFrames(i));
    pressure = p(i,1:StudyFrames(i));
    
    if StudyFrames(i) == 30
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
    volume = v(i,1:StudyFrames(i));
    pressure = p(i,1:StudyFrames(i));
    
    % Get diastasis point (minimum in LVP)
    ds = find(pressure == min(pressure));
    
    if StudyFrames(i) == 30
        h1 = plot(volume(ds:end), pressure(ds:end), 'm-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
    else
        h2 = plot(volume(ds:end), pressure(ds:end), 'k-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
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
    volume = v(i,1:StudyFrames(i));
    pressure = p(i,1:StudyFrames(i));
    
    % Get diastasis point (minimum in LVP)
    ds = find(pressure == min(pressure));
    
    if StudyFrames(i) == 30
        h1 = plot(volume([ds,end]), pressure([ds,end]), 'm-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');
    else
        h2 = plot(volume([ds,end]), pressure([ds,end]), 'k-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end
    hold on
    
end

xlabel('Volume (mL)', 'FontSize', 12)
ylabel('Pressure (kPa)', 'FontSize', 12)
legend([h1, h2], {'Study 1', 'Study 2'}, 'Location', 'Best')
saveas(FH3, sprintf('%s/PV-Diastole-DS-ED-Study%s.png', parentDir, descriptor));

close all