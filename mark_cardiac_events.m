function [ diastasis_index, endIVR_index, ED_index ] = mark_cardiac_events( representative_LV_pressure, deriv1_LV, deriv2_LV )

% End IVR - point of inflexion and negative first derivative (slope greater
% than -150 to ignore noise and flat regions)
modified_deriv2LV = zeros(1, length(deriv2_LV));
for i = 1:length(deriv2_LV)
    if ((deriv1_LV(i)<=-150) && (deriv2_LV(i)>0))
        modified_deriv2LV(i) = deriv2_LV(i);
    else
        modified_deriv2LV(i) = min(deriv2_LV);
    end
end 

endIVR_index = find(max(modified_deriv2LV) == modified_deriv2LV); %%% RM - I don't understand this.... I thought that to find the point of inflection, you find where the second derivative is equal to zero, not the maximum.....
endIVR_index = endIVR_index(1); %% RM - are there multiple? Why only take the first one??

% End diastole - R peak of ECG so last point of representative trace;
ED_index = length(representative_LV_pressure);

% Diastasis - rate of change of pressure is 0 and between Pmin and ED
min_pressure_index = (find(min(representative_LV_pressure)==representative_LV_pressure));
DS_deriv2LV = zeros(1,length(deriv2_LV));
for i = 1:length(deriv2_LV)
    if ((i>min_pressure_index) && (i<ED_index))
        DS_deriv2LV(i) = deriv2_LV(i);
    else
        DS_deriv2LV(i) = max(deriv2_LV);
    end
end

% Diastasis after first change of sign in second derivative from positive to
% negative - at minimum second derivative is positive
for i = min_pressure_index:ED_index
    if DS_deriv2LV(i)<0
        diastasis_index=i;
        break
    end
end
    
diastasis_index = diastasis_index(1);

end

