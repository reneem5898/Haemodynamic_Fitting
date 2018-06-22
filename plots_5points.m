function plots_5points(timepoints_MRI, mri_time, LVP, ed2eivc, eivc2es, es2eivr, eivr2ds, ds2ed, V, directory)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates plots of resulting interpolated pressures
%
% Inputs: 1) timepoints_MRI: vector of cardiac timepoints from MRI
%         2) mri_time: vector of time of one cardiac cycle
%         3) LVP: final averaged left ventricular pressure trace
%         4) ed2eivc: LV presure between ED and eIVC
%         5) eivc2es: LV pressure between eIVC and ES
%         6) es2eivr: LV pressure between ES and eIVR
%         7) eivr2ds: LV pressure between eIVR and DS
%         8) ds2ed: LV pressure between DS and ED
%         9) V: volumes (from CIM model)
%         10) directory: directory where to save plots
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
eivcMRI = timepoints_MRI(2);
esMRI = timepoints_MRI(3);
eivrMRI = timepoints_MRI(4);
dsMRI = timepoints_MRI(5);


%%%%%%%%%%%%%%%% Beat-averaged left ventricular pressure %%%%%%%%%%%%%%%%%%
FH8 = figure;
% Plot ED to eIVC
t = repmat(mri_time(edMRI:eivcMRI),size(ed2eivc,1),1);
scatter(t(:), ed2eivc(:), 'g*')
hold on
% Plot eIVC to ES
t = repmat(mri_time((eivcMRI+1):esMRI),size(eivc2es,1),1);
scatter(t(:), eivc2es(:), 'r*')
% Plot ES to eIVR
t = repmat(mri_time((esMRI+1):eivrMRI),size(es2eivr,1),1);
scatter(t(:), es2eivr(:), 'm*')
% Plot eIVR to DS
t = repmat(mri_time((eivrMRI+1):dsMRI-1),size(eivr2ds,1),1);
scatter(t(:), eivr2ds(:), 'c*')
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
% Plot ED to eIVC
t = repmat(V(edMRI:eivcMRI)',size(ed2eivc,1),1);
scatter(t(:), ed2eivc(:), 'g*')
hold on
% Plot eIVC to ES
t = repmat(V((eivcMRI+1):esMRI)',size(eivc2es,1),1);
scatter(t(:), eivc2es(:), 'r*')
% Plot ES to eIVR
t = repmat(V((esMRI+1):eivrMRI)',size(es2eivr,1),1);
scatter(t(:), es2eivr(:), 'm*')
% Plot eIVR to DS
t = repmat(V((eivrMRI+1):dsMRI-1)',size(eivr2ds,1),1);
scatter(t(:), eivr2ds(:), 'c*')
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
scatter(mri_time(eivrMRI), V(eivrMRI), 50, 'ro', 'filled')
scatter(mri_time(eivcMRI), V(eivcMRI), 50, 'ro', 'filled')
xlabel('Time (s)', 'FontSize', 16)
ylabel('Volume (mL)', 'FontSize', 16)
set(gca, 'FontSize', 12)
saveas(FH10, sprintf('%s/volume-mri-cardiac-events.png', directory));