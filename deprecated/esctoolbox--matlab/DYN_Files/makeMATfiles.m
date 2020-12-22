% makeMATfiles:
% This utility script loads the Excel files produced by the Arbin cell 
% tester and converts them to MATLAB ".mat" files for easier access.

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.


setupDynData; % get list of files to be processed
% Column headers to look for and convert to ".mat" file
headers = {'Test_Time(s)','Step_Index','Current(A)','Voltage(V)',...
           'Charge_Capacity(Ah)','Discharge_Capacity(Ah)'};
% Corresponding MATLAB structure field names to use         
finFields = {'time','step','current','voltage','chgAh','disAh'};
% Field names to use for the four different testing scripts
stepFields = {'script1','script2'};
% Set skipDone to "1" to skip processing if ".mat" file already exists, or
% to "0" to reprocess file
skipDone = 0;

for indID = 1:length(cellIDs),      % loop over all cell types
  for indTemps = 1:length(temps),   % loop over all temperatures
    if indTemps > length(mags{indID}), break, end % skip if no data
    theMag = mags{indID}(indTemps); % relative C-rate of data file
    if theMag < 0, continue, end    % skip if no data
    data = [];                      % clear data structure and start fresh
    if temps(indTemps) < 0, % Use this filename for negative temps
      DYNPrefix = sprintf('%s_DYN/%s_DYN_%02d_N%02d',...
        cellIDs{indID},cellIDs{indID},theMag,abs(temps(indTemps)));
    else                    % Use this filename for positive temps
      DYNPrefix = sprintf('%s_DYN/%s_DYN_%02d_P%02d',...
        cellIDs{indID},cellIDs{indID},theMag,temps(indTemps));
    end
    OUTFile = sprintf('%s.mat',DYNPrefix); % output filename, incl. path
    if exist(OUTFile,'file') && skipDone, 
      fprintf('Skipping %s: already done\n',OUTFile); 
      continue
    end
    DYNFile1 = sprintf('%s_S1.xlsx',DYNPrefix); % input file, script 1
    DYNFile2 = sprintf('%s_S2.xlsx',DYNPrefix); % input file, script 2
    if ~exist(DYNFile1,'file'), 
      fprintf('Skipping %s: Missing source data\n',DYNFile1); 
      continue
    end
    if ~exist(DYNFile2,'file'), 
      fprintf('Skipping %s: Missing source data\n',DYNFile2); 
      continue
    end
      
    for theScript = 1:2, % process both scripts
      scriptData = [];   % clear data structure for script
      for theField = 1:length(finFields), % initialize output fields
        scriptData.(finFields{theField}) = [];
      end

      DYNFile = sprintf('%s_S%d.xlsx',DYNPrefix,theScript); 
      [~,sheets] = xlsfinfo(DYNFile); % get worksheet names
      fprintf('Reading %s\n',DYNFile); 
      for theSheet = 1:length(sheets), % process all worksheets
        if strcmp(sheets{theSheet},'Info'), continue; end % except "Info"
        fprintf('  Processing sheet %s\n',sheets{theSheet});
        [num,txt,raw] = xlsread(DYNFile,sheets{theSheet}); % read data
        for theHead = 1:length(headers), % scan for desired data
          ind = strcmp(txt(1,:),headers{theHead}); % it's in this column
          scriptData.(finFields{theHead}) = ... % store data
            [scriptData.(finFields{theHead}); num(:,ind == 1)];
        end
      end
      DYNData.(stepFields{theScript}) = scriptData; % save data
    end % for theScript

    % Do some minor processing: 1-s interpolation
    ind = find(diff(DYNData.script1.time)<=0); % delete data with duplicate
    DYNData.script1.time(ind+1)=[];            % time stamp so "interp1" 
    DYNData.script1.voltage(ind+1)=[];         % won't get confused
    DYNData.script1.current(ind+1)=[]; 
    DYNData.script1.chgAh(ind+1)=[];
    DYNData.script1.disAh(ind+1)=[];
    DYNData.script1.step(ind+1)=[];

    % Interpolate raw measured data in 1-s increments
    ind = find(DYNData.script1.step == 2,1); % keep step 2, and last 300s 
    t1=DYNData.script1.time(ind) - 300; t2=DYNData.script1.time(end); %of
    DYNData.script1.current = -interp1(DYNData.script1.time,... % step 1
      DYNData.script1.current,t1:t2,'linear'); % 1-s i(t)
    DYNData.script1.voltage = interp1(DYNData.script1.time,...
      DYNData.script1.voltage,t1:t2,'linear'); % 1-s v(t)
    DYNData.script1.chgAh = interp1(DYNData.script1.time,...
      DYNData.script1.chgAh,t1:t2,'linear'); % 1-s v(t)
    DYNData.script1.disAh = interp1(DYNData.script1.time,...
      DYNData.script1.disAh,t1:t2,'linear'); % 1-s v(t)
    DYNData.script1.step = interp1(DYNData.script1.time,...
      DYNData.script1.step,t1:t2,'nearest'); % 1-s v(t)
    DYNData.script1.time = t1:t2;

    % Same for script-2 data
    ind = find(diff(DYNData.script2.time)<=0); 
    DYNData.script2.time(ind+1)=[];    
    DYNData.script2.voltage(ind+1)=[]; 
    DYNData.script2.current(ind+1)=[]; 
    DYNData.script2.chgAh(ind+1)=[];
    DYNData.script2.disAh(ind+1)=[];
    DYNData.script2.step(ind+1)=[];
    t1=DYNData.script2.time(1); t2=DYNData.script2.time(end);    
    DYNData.script2.current = -interp1(DYNData.script2.time,...
      DYNData.script2.current,t1:t2,'linear'); % 1-s i(t)
    DYNData.script2.voltage = interp1(DYNData.script2.time,...
      DYNData.script2.voltage,t1:t2,'linear'); % 1-s v(t)
    DYNData.script2.chgAh = interp1(DYNData.script2.time,...
      DYNData.script2.chgAh,t1:t2,'linear'); % 1-s v(t)
    DYNData.script2.disAh = interp1(DYNData.script2.time,...
      DYNData.script2.disAh,t1:t2,'linear'); % 1-s v(t)
    DYNData.script2.step = interp1(DYNData.script2.time,...
      DYNData.script2.step,t1:t2,'nearest'); % 1-s v(t)
    DYNData.script2.time = t1:t2;

    % Do some minor processing: split script2 into two parts
    % look for longest period of charging; split just before this
    % flagI is 0 if discharging or "step" if charging
    flagI = (DYNData.script2.current < 0).*DYNData.script2.step; 
    starts = find([1 diff(flagI)]); % find starting point of each run
    runs = diff(find(diff(flagI))); % find length of all runs
    splitInd = starts(find(runs == max(runs))+1);

    DYNData.script3.time = DYNData.script2.time(splitInd:end);
    DYNData.script2.time(splitInd:end) = [];
    DYNData.script3.voltage = DYNData.script2.voltage(splitInd:end);
    DYNData.script2.voltage(splitInd:end) = [];
    DYNData.script3.current = DYNData.script2.current(splitInd:end);
    DYNData.script2.current(splitInd:end) = [];
    DYNData.script3.chgAh = DYNData.script2.chgAh(splitInd:end) - ...
                            DYNData.script2.chgAh(splitInd);    
    DYNData.script2.chgAh(splitInd:end) = [];
    DYNData.script3.disAh = DYNData.script2.disAh(splitInd:end) - ...
                            DYNData.script2.disAh(splitInd);
    DYNData.script2.disAh(splitInd:end) = [];
    DYNData.script3.step = DYNData.script2.step(splitInd:end);
    DYNData.script2.step(splitInd:end) = [];

    save(OUTFile,'DYNData'); % save processed data
  end % for indTemps
end % for indID