% --------------------------------------------------------------------
% script runProcessDynamic
%
% RUNPROCESSDYNAMIC reads data files corresponding to dynamic cell 
% tests, executes PROCESSDYNAMIC, and then saves the resulting 
% ESC model.  It relies on SETUPDYNDATA to provide a list of data files
% to be processed.

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

close all
setupDynData; % get list of files to be processed
numpoles = 1; % number of resistor--capacitor pairs in final model

for indID = 1:length(cellIDs), % process each cell type
  cellID = cellIDs{indID};     % get cell identifier
  
  % Read model OCV file, previously computed by runProcessOCV
  modelFile = sprintf('../OCV_Files/%smodel-ocv.mat',cellID);
  if ~exist(modelFile,'file'),
    error(['File "%s" not found.\n' ...
      'Please change folders so that "%s" points to a valid model '...
      'file and re-run runProcessDynamic.'],modelFile,modelFile); %#ok<SPERR>
  end
  load(modelFile);
  
  % Read MAT raw data files
  data = zeros([0 length(mags{indID} > 0)]); dataInd = 0;
  for indTemps = 1:length(mags{indID}), % read all temperatures
    theMag = mags{indID}(indTemps);     % max C-rate in data file * 10
    if theMag < 0,                      % omit these data files
      continue 
    else                                % store this data in "data"
      dataInd = dataInd + 1;
    end
    if temps(indTemps) < 0, % if temperature is negative, then load this
      DYNPrefix = sprintf('%s_DYN/%s_DYN_%02d_N%02d',... % data file
        cellID,cellID,theMag,abs(temps(indTemps)));
    else                    % if temperature is positive, then load this
      DYNPrefix = sprintf('%s_DYN/%s_DYN_%02d_P%02d',... % data file
        cellID,cellID,theMag,temps(indTemps));
    end
    inFile = sprintf('%s.mat',DYNPrefix);
    if ~exist(inFile,'file'),
      error(['File "%s" not found.\n' ...
        'Please change folders so that "%s" points to a valid data '...
        'file and re-run runProcessDynamic.'],inFile,inFile); %#ok<SPERR>
    end
    fprintf('Loading %s\n',inFile); load(inFile);        
    data(dataInd).temp    = temps(indTemps); % store temperature
    data(dataInd).script1 = DYNData.script1; % store data from each of the
    data(dataInd).script2 = DYNData.script2; % three scripts
    data(dataInd).script3 = DYNData.script3;
  end
  
  model = processDynamic(data,model,numpoles,1); % does the "heavy lifting"
  modelFile = sprintf('%smodel.mat',cellID); % save optimized model in this
  save(modelFile,'model');                   % file
  
  % Plot model-match voltage results at 25 degC, plus RMS voltage-estimation 
  % error between 5% and 95% cell state of charge
  figure(10+indID);
  indTemps = find(temps == 25);
  [vk,rck,hk,zk,sik,OCV] = simCell(data(indTemps).script1.current,...
    temps(indTemps),1,model,1,zeros(numpoles,1),0);
  tk = (1:length(data(indTemps).script1.current))-1;
  plot(tk,data(indTemps).script1.voltage,tk,vk);
  verr = data(indTemps).script1.voltage - vk';
  v1 = OCVfromSOCtemp(0.95,temps(indTemps),model);
  v2 = OCVfromSOCtemp(0.05,temps(indTemps),model);
  N1 = find(data(indTemps).script1.voltage<v1,1,'first'); 
  N2 = find(data(indTemps).script1.voltage<v2,1,'first');
  if isempty(N1), N1=1; end; if isempty(N2), N2=length(verr); end
  rmserr=sqrt(mean(verr(N1:N2).^2));
  fprintf('RMS error of simCell @ 25 degC = %0.2f (mv)\n',rmserr*1000);
end