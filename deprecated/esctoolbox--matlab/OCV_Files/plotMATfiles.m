% This utility script loads the MATLAB ".mat" files produced by 
% makeMATfiles.m and plots the voltages recorded for the four scripts

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

clear all; close all

cellIDs = {'A123','ATL','E1','E2','P14','SAM'}; % Identifiers for each cell
order = [-25 -15 -5 5 15 25 35 45];             % Temperatures for each 
% Field names to use for the four different testing scripts
stepFields = {'script1','script2','script3','script4'};

for theID = 1:length(cellIDs),    % loop over all cells
  data = [];                      % clear data structure and start fresh
  for theFile = 1:length(order),  % loop over all temperatures
    if order(theFile) < 0,        % if temperature is negative, then
      OCVPrefix = sprintf('%s_OCV/%s_OCV_N%02d',... % load this file
        cellIDs{theID},cellIDs{theID},abs(order(theFile)));
    else                          % if temperature is positive, then
      OCVPrefix = sprintf('%s_OCV/%s_OCV_P%02d',... % load this file
        cellIDs{theID},cellIDs{theID},order(theFile));
    end
    inFile = sprintf('%s.mat',OCVPrefix); % full filename, including path
    if ~exist(inFile,'file'),
      error(['File "%s" not found in current folder.\n' ...
        'Please change folders so that "%s" is in the current '...
        'folder and re-run plotMATfiles.'],inFile,inFile); %#ok<SPERR>
    end    
    load(inFile);                 % load the file

    figure                        % make new figure
    for theScript = 1:4,          % plot each script's data
      subplot(2,2,theScript);
      data = OCVData.(stepFields{theScript});
      plot(data.time,data.voltage);
      title(sprintf('%s - step %d',OCVPrefix,theScript),'interpreter','none');
    end
  end
end