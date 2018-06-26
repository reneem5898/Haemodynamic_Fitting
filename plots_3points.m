function plots_3points(timepoints_MRI, mri_time, LVP, ed2es, es2ds, ds2ed, V, directory)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates plots of resulting interpolated pressures
%
% Inputs: 1) timepoints_MRI: vector of cardiac timepoints from MRI
%         2) mri_time: vector of time of one cardiac cycle
%         3) LVP: final averaged left ventricular pressure trace
%         4) ed2es: LV presure between ED and ES
%         5) es2ds: LV pressure between ES and DS
%         6) ds2ed: LV pressure between DS and ED
%         7) V: volumes (from CIM model)
%         8) directory: directory where to save plots
%
% Plots: 1) beat averaged left ventricular pressure
%        2) pressure volume loop
%        3) volume curve with highlighted points
%
% Written by: Renee Miller (renee.miller@auckland.ac.nz)
% Date: 22 June 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MRI cardiac time points
edMRI = timepoints_MRI(1);
esMRI = timepoints_MRI(2);
dsMRI = timepoints_MRI(3);


%%%%%%%%%%%%%%%% Beat-averaged left ventricular pressure %%%%%%%%%%%%%%%%%%
FH8 = figure;
% Plot ED to ES
t = repmat(mri_time(edMRI:esMRI),size(ed2es,1),1);
scatter(t(:), ed2es(:), 'g*')
hold on
% Plot ES to DS
t = repmat(mri_time((esMRI+1):(dsMRI-1)),size(es2ds,1),1);
scatter(t(:), es2ds(:), 'm*')
% Plot DS to ED
if edMRI == 1
    t = repmat(mri_time(dsMRI:end),size(ds2ed,1),1);
else
    t = repmat([mri_time(dsMRI:end),mri_time(1:(edMRI-1))],size(ds2ed,1),1);
end
scatter(t(:), ds2ed(:), 'b*')
% Plot average
plot(mri_time, LVP, 'k*-')
xlabel('Time (s)', 'FontSize', 16)
ylabel('LV Pressure (kPa)', 'FontSize', 16)
title('Beat-averaged LV Pressure', 'FontSize', 16)
set(gca, 'FontSize', 12)
saveas(FH8, sprintf('%s/beat-averaged-lv-pressure.png', directory));


%%%%%%%%%%%%%%%%%%% Pressure-volume loop %%%%%%%%%%%%%%%%%%
FH9 = figure;
% Plot ED to ES
t = repmat(V(edMRI:esMRI)',size(ed2es,1),1);
scatter(t(:), ed2es(:), 'g*')
hold on
% Plot ES to eIVR
t = repmat(V((esMRI+1):(dsMRI-1))',size(es2ds,1),1);
scatter(t(:), es2ds(:), 'm*')
% Plot DS to ED
if edMRI == 1
    t = repmat(V(dsMRI:end)',size(ds2ed,1),1);
else
    t = repmat([V(dsMRI:end),V(1:(edMRI-1))]',size(ds2ed,1),1);
end
scatter(t(:), ds2ed(:), 'b*')
plot(V, LVP, 'k*-')
xlabel('Volume (mL)', 'FontSize', 16)
ylabel('LV Pressure (kPa)', 'FontSize', 16)
title('Pressure-Volume Loop', 'FontSize', 16)
set(gca, 'FontSize', 12)
saveas(FH9, sprintf('%s/pv-loop.png', directory));


%%%%%%%%%%%%%%%%%%% Volume with cardiac events %%%%%%%%%%%%%%%%%%
FH10 = figure;
plot(mri_time, V, 'ksq-', 'LineWidth', 2, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
hold on
scatter(mri_time(dsMRI), V(dsMRI), 50, 'ro', 'filled')
scatter(mri_time(edMRI), V(edMRI), 50, 'ro', 'filled')
scatter(mri_time(esMRI), V(esMRI), 50, 'ro', 'filled')
xlabel('Time (s)', 'FontSize', 16)
ylabel('Volume (mL)', 'FontSize', 16)
set(gca, 'FontSize', 12)
saveas(FH10, sprintf('%s/volume-mri-cardiac-events.png', directory));