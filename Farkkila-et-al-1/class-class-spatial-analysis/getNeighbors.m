function [lengthList, neighbour_columns] =getNeighbors(xydata, max_dist)
%GETNEIGHBORS creates new columns with the neighbor IDs in similar format
%to histoCAT output.
%
%  xydata - matrix of datapoints, each row represents one cell
%  max_dist -  in pixels for largest search
% Blame to Julia Casado

%% Set default max_dist if not included in function call
if nargin<2
    max_dist = 28;
end
radius = max_dist;

%% Find neighbors wihin the defined radius
cellarray_neighbourids = rangesearch(xydata, xydata, radius);

%% Convert neighbors cellarray to table
lengthList = cellfun('size',cellarray_neighbourids,2);
maxLength = max(lengthList);
neighbour_CellId = zeros(length(cellarray_neighbourids),maxLength);
for i=1:length(cellarray_neighbourids)
    neighbour_CellId(i,1:lengthList(i)) = cellarray_neighbourids{i,1};
end

%% Prettify table to match HistoCAT output
neighbour_names = strcat({strcat('neighbour_5_CellId')}, int2str((1:maxLength).')).';
neighbour_names_include = strrep(neighbour_names,' ','');

neighbour_columns = array2table(neighbour_CellId,'VariableNames',neighbour_names_include);

end