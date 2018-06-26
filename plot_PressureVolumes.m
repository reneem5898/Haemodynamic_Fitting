function plot_PressureVolumes(p, v, parentDir, descriptor, omit, TimePoints, StudyFrames)

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
    
    if strcmp(TimePoints{i},'Base') && ~any(omit == i)
        h1 = plot(volume, pressure, 'r-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
        hold on
    elseif strcmp(TimePoints{i},'M1') && ~any(omit == i)
        h2 = plot(volume, pressure, 'b-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
        hold on
    elseif strcmp(TimePoints{i},'M2') && ~any(omit == i)
        h3 = plot(volume, pressure, 'g-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
        hold on
    end
end
xlabel('Volume (mL)', 'FontSize', 12)
ylabel('Pressure (kPa)', 'FontSize', 12)
legend([h1, h2, h3], {'Bx', 'M1', 'M2'}, 'Location', 'Best')
saveas(FH, sprintf('%s/PV-Loops%s.png', parentDir, descriptor));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot all diastolic portions of the PV loops together
FH2 = figure;
for i = 1:numCases
    
    % Get pressure and volume for current study
    volume = v(i,1:StudyFrames(i));
    pressure = p(i,1:StudyFrames(i));
    
    % Get diastasis point (minimum in LVP)
    ds = find(pressure == min(pressure));
    
    % Plot pressure between diastasis and end diastole
    if strcmp(TimePoints{i},'Base') && ~any(omit == i)
        h1 = plot(volume(ds:end), pressure(ds:end), 'r-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
        hold on
    elseif strcmp(TimePoints{i},'M1') && ~any(omit == i)
        h2 = plot(volume(ds:end), pressure(ds:end), 'b-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
        hold on
    elseif strcmp(TimePoints{i},'M2') && ~any(omit == i)
        h3 = plot(volume(ds:end), pressure(ds:end), 'g-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
        hold on
    end
end
xlabel('Volume (mL)', 'FontSize', 12)
ylabel('Pressure (kPa)', 'FontSize', 12)
legend([h1, h2, h3], {'Bx', 'M1', 'M2'}, 'Location', 'Best')
saveas(FH2, sprintf('%s/PV-Diastole%s.png', parentDir, descriptor));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot DS and ED
FH3 = figure;
for i = 1:numCases
    
    % Get pressure and volume for current study
    volume = v(i,1:StudyFrames(i));
    pressure = p(i,1:StudyFrames(i));
    
    % Get diastasis point (minimum in LVP)
    ds = find(pressure == min(pressure));
    
    % Plot pressure at diastasis and end diastole only
    if strcmp(TimePoints{i},'Base') && ~any(omit == i)
        h1 = plot(volume([ds,end]), pressure([ds,end]), 'r-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
        hold on
    elseif strcmp(TimePoints{i},'M1') && ~any(omit == i)
        h2 = plot(volume([ds,end]), pressure([ds,end]), 'b-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
        hold on
    elseif strcmp(TimePoints{i},'M2') && ~any(omit == i)
        h3 = plot(volume([ds,end]), pressure([ds,end]), 'g-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
        hold on
    end
end

xlabel('Volume (mL)', 'FontSize', 12)
ylabel('Pressure (kPa)', 'FontSize', 12)
legend([h1, h2, h3], {'Bx', 'M1', 'M2'}, 'Location', 'Best')
saveas(FH3, sprintf('%s/PV-Diastole-DS-ED%s.png', parentDir, descriptor));

close all