% This utility script loads the MATLAB ".mat" files produced by 
% makeMATfiles.m and plots the voltages recorded for the three scripts

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

close all

% Set up cell identifiers and multipliers at different temperatures
setupDynData
% Field names to use for the three different testing scripts
stepFields = {'script1','script2','script3'}; % 

for indID = 1:length(cellIDs), % loop over all cells
  cellID = cellIDs{indID};

  % Read MAT raw data files
  data = zeros([0 length(mags{indID} > 0)]); dataInd = 0;
  for indTemps = 1:length(mags{indID}), % loop over all temperatures
    theMag = mags{indID}(indTemps);     % if data does not exist, skip
    if theMag < 0, 
      continue 
    else                                % store data at this index
      dataInd = dataInd + 1;
    end
    if temps(indTemps) < 0,             % if temperature is negative, load
      DYNPrefix = sprintf('%s_DYN/%s_DYN_%02d_N%02d',... % this file
        cellID,cellID,theMag,abs(temps(indTemps)));
    else                                % if temperature is positive, load
      DYNPrefix = sprintf('%s_DYN/%s_DYN_%02d_P%02d',... % this file
        cellID,cellID,theMag,temps(indTemps));
    end
    inFile = sprintf('%s.mat',DYNPrefix); % full filename with path
    if ~exist(inFile,'file'),
      error(['File "%s" not found in current folder.\n' ...
        'Please change folders so that "%s" is in the current '...
        'folder and re-run plotMATfiles.'],inFile,inFile); %#ok<SPERR>
    end    
    fprintf('Loading %s\n',inFile); load(inFile);        

    figure
    for theScript = 1:3, % loop over all scripts
      subplot(2,2,theScript);
      data = DYNData.(stepFields{theScript});
      t = (data.time - data.time(1))/3600;
      plot(t,data.voltage); % plot the voltage for this script
      title(sprintf('%s - step %d',DYNPrefix,theScript),'interpreter','none');
    end
  end
end