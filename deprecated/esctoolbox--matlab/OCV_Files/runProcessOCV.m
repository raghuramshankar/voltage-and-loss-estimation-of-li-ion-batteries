% script runProcessOCV.m
%   Loads data from OCV lab tests done for several cells at several 
%   different temperatures, then calls processOCV.m to create the OCV
%   relationship, then saves the model to a model file.

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

clear all
cellIDs = {'A123','ATL','E1','E2','P14','SAM'}; % Identifiers for each cell
% data files for each cell available at these temperatures
temps = {[-25 -15 -5 5 15 25 35 45], ... % A123
         [-25 -15 -5 5 15 25 35 45], ... % ATL
         [-25 -15 -5 5 15 25 35 45], ... % E1
         [-25 -15 -5 5 15 25 35 45], ... % E2
         [-25 -15 -5 5 15 25 35 45], ... % P14
         [-25 -15 -5 5 15 25 35 45]};    % SAM
% minimum and maximum voltages for each cell, used for plotting results       
minV = [2.00 2.00 2.95 2.50 2.50 2.50];
maxV = [3.75 3.75 4.25 4.25 4.25 4.25];

% --------------------------------------------------------------------
% Load raw data from cell tests, then process
% --------------------------------------------------------------------
for theID = 1:length(cellIDs), % loop over all cells
  dirname = cellIDs{theID}; cellID = dirname;
  ind = find(dirname == '_'); % if there is a "_", delete it
  if ~isempty(ind), dirname = dirname(1:ind-1); end
  OCVDir = sprintf('%s_OCV',dirname); % folder in which to find data
  if ~exist(OCVDir,'dir'),
    error(['Folder "%s" not found in current folder.\n' ...
      'Please change folders so that "%s" is in the current '...
      'folder and re-run runProcessOCV.'],OCVDir,OCVDir); %#ok<SPERR>
  end
  
  filetemps = temps{theID}(:);  % data exists at these temperatures
  numtemps = length(filetemps); % number of data sets
  data = zeros([0 numtemps]);   % initialize data to zero

  for k = 1:numtemps,           % load the data files into the "data" var
    if filetemps(k) < 0,        % if temperature is negative, then
      filename = sprintf('%s/%s_OCV_N%02d.mat',... % look for this file
        OCVDir,cellID,abs(filetemps(k)));
    else                        % if temperature is positive, then
      filename = sprintf('%s/%s_OCV_P%02d.mat',... % look for this file
        OCVDir,cellID,filetemps(k));
    end
    load(filename);             % load OCV data file
    data(k).temp = filetemps(k);       % save temperature of test
    data(k).script1 = OCVData.script1; % save the four scripts
    data(k).script2 = OCVData.script2;
    data(k).script3 = OCVData.script3;
    data(k).script4 = OCVData.script4;
  end

  % then, call "processOCV" to do the actual data processing
  model = processOCV(data,cellID,minV(theID),maxV(theID),1);
  save(sprintf('%smodel-ocv.mat',cellID),'model'); % save model file
end