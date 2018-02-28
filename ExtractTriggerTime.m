function [RR_mean, TS_mean, TS_std, no_frames] = ExtractTriggerTime(image_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is designed to extract trigger time from all dicom images
%
% Originally written by: Jenny Wang
% Modified by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 13 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get list of files in the image directory
AllFiles = dir(image_dir);
files = {AllFiles.name};
files(ismember(files,{'.','..'})) = [];

% Get number of frames from header of one image file
info = dicominfo(strcat([image_dir,'\',files{1}]));
no_frames = info.CardiacNumberOfImages;

no_slices = length(files)/no_frames; % Calculate number of slices
k = 1; % Counter
TT_all_frames = zeros(no_frames, no_slices); % Variable to save trigger times
TS = zeros(no_slices,1); % Variable to save temporal spacing
sliceLocations = zeros(length(files),1);


%% Get list of slice locations - not always collected in order
for i = 1:no_slices % Loop through number of slices
    for j = 1:no_frames % Loop through number of frames
        infor = dicominfo(strcat([image_dir,'\',files{k}])); % Get info from header
        sliceLocations(k) = infor.SliceLocation;
        k = k+1; % Increment image number
    end
end

% Get unique list of slice locations in order
sliceLocationsUnique = unique(sliceLocations);

k = 1; % counter
for i = 1:no_slices
    
    for j = 1:no_frames % Loop through number of frames
        infor = dicominfo(strcat([image_dir,'\',files{k}])); % Get info from header
        % Get index of slice location - not always collected in order      
        TT_all_frames(j,sliceLocationsUnique == infor.SliceLocation) = infor.TriggerTime; % Save all trigger times
        k = k+1; % Increment image number
    end
    
    %% Sort the trigger time to see whether there is an image out of place (?)
    [TT_sort, index] = sortrows(TT_all_frames(:,i));
    frame_number = linspace(1, no_frames, no_frames);
    frame_diff = index - frame_number';
    if any(frame_diff~=0)
        fprintf('*** Warning: Slice %d requires reordering of the frames ......\n',i);
        TT_all_frames(:,i) = TT_sort;
    end
    
    %% Calculate the temporal spacing
    TS_all_frames = diff(TT_sort);
    if diff(TS_all_frames) < 1e-4
        TS(i) = TS_all_frames(1);
        fprintf('*** Slice %d has consistent temporal spacing ......\n',i);
    else
        fprintf('*** Warning: Slice %d does not have consistent temporal spacing ......\n',i); % RM - Then what happens??
    end
end

%% Calculate the mean temporal spacing
TS_mean = mean(TS);
TS_std = std(TS);
RR_mean = TS_mean*(no_frames-1);

fprintf('+++++ The mean temporal spacing is %f .....\n', TS_mean);
fprintf('+++++ The standard deviation for temporal spacing is %f .....\n', TS_std);
fprintf('+++++ The mean R-R interval is %f .........\n', RR_mean);

return
