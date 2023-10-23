% getEnergy – code returns delta power in a given band and power spectra.
% also saves spectral responses of individual protocols

% Inputs: subjectName ,
%         expDate,
%         protocolName,
%         freqRange
%         tBL,tST - time range for baseline, and stimulus 
%         electrodeNum - scalar electrode number or vector of elecs
% outputs: logDelPwrBand: the gamma power change from baseline for that electrode. log10 
%          meanPwrST:    power in time range of ST
%          meanPwrBLall: power in time range of BL
%          freqVals:     frequency axis
%          parameterCombinations: stimulus conditions for the protocol         
% STK 240822 

function [delPwrBand,meanPwrST,meanPwrBL,bandPwrST,bandPwrBL,parameterCombinations,freqVals,delPwrBandF,meanPwrSTF,meanPwrBLF,bandPwrSTF,bandPwrBLF]= getEnergy(subjectName,expDate,protocolName,fBandGamma,fBad,tBL,tST,electrodeNum,numTapers)

gridType            = 'Microelectrode';
folderSourceString  =pwd; % folder where 'data' directory is kept with extracted protocols
if ~exist('fBandGamma','var'),             fBandGamma = [30 80]; end
if ~exist('fBad','var'),                   fBad = [48:52, 98:102, 148:152]; end
if ~exist('tBL','var'),                    tBL = [-0.25 0];  end
if ~exist('tST','var'),                    tST = [0.25 0.5]; end
if ~exist('electrodeNum','var'),           electrodeNum = 2; end
if ~exist('numTapers','var'),              numTapers = 1;    end  
% used for Multitaper by Chronux
modelDataDir = fullfile('savedData','spectra');% fullfile('modelData','derivatives','spectra');  % where saved power response files are kept
if ~exist(modelDataDir,'dir'), mkdir(modelDataDir);   end

% get the responses for this particular experiment:
[powerST,powerBL,freqVals,electrodes,parameterCombinations] = getResponses(gridType,folderSourceString,subjectName,expDate,protocolName,tBL,tST,numTapers,modelDataDir);

% select the electrodes:
elPowerST = cell(length(electrodeNum),size(powerST,2),size(powerST,3),size(powerST,4),size(powerST,5),size(powerST,6));
elPowerBL = cell(length(electrodeNum),size(powerBL,2),size(powerBL,3),size(powerBL,4),size(powerST,5),size(powerST,6));
[~,eind]  = intersect(electrodes,electrodeNum);
elPowerST(:,:,:,:,:,:) = powerST(eind,:,:,:,:,:,:); 
elPowerBL(:,:,:,:,:,:) = powerBL(eind,:,:,:,:,:,:);

% get band power and trial averaged spectra
[delPwrBand,meanPwrST,meanPwrBL,bandPwrST,bandPwrBL] = getPowerChanges(elPowerST,elPowerBL,freqVals,{fBandGamma},fBad);
  
% if required for different folds?
[~,numTrials,~]  = cellfun(@size, elPowerBL(1,:,:,:,:,:));
ind = find(numTrials==1);
for nt=1:length(ind)
    elPowerBL{ind(nt)}=repmat(elPowerBL{ind(nt)},[1 2]);   % repeat the same trial for folds
    elPowerST{ind(nt)}=repmat(elPowerST{ind(nt)},[1 2]);
    [~,numTrials,~]  = cellfun(@size, elPowerBL(1,:,:,:,:,:));
end
numTrialsC = cellfun(@(x) cvpartition(x,'k',2),num2cell(numTrials),'un',0);
numTrialsFold= repmat(numTrialsC,size(elPowerST,1),1,1);
for f = 1:2    % 2 folds made
    elPowerSTF = cellfun(@(x,y) x(:,y.training(f),:), elPowerST,numTrialsFold,'un',0);
    elPowerBLF = cellfun(@(x,y) x(:,y.training(f),:), elPowerBL,numTrialsFold,'un',0);
    [delPwrBandF{f},meanPwrSTF{f},meanPwrBLF{f},bandPwrSTF{f},bandPwrBLF{f}] = getPowerChanges(elPowerSTF,elPowerBLF,freqVals,{fBandGamma},fBad);
end  
end

% % 
function [pwrST,pwrBL,freqVals,electrodeList,parameterCombinations]= getResponses(gridType,folderSourceString,subjectName,expDate,protocolName,tBL,tST,numTapers,saveDir)

if ~exist('numTapers','var')||isempty(numTapers), numTapers = 3; end

% 1. load good electrodes for this monkey
[~,~,LFPElectrodeList,EcogElectrodeList,~] = getRFdetails(subjectName,'savedData'); % savedData dir has 'RFDetails.mat' file
electrodeList=[LFPElectrodeList{1}; EcogElectrodeList{1}]; % unwrap them from cell.

% 2. get responses for protocol. returns in terms of elec x size x sf x con x ori x tf
if ~isempty(electrodeList)
      [pwrST,pwrBL,freqVals,parameterCombinations]  = getPowerMT(folderSourceString,gridType,subjectName,expDate,protocolName,electrodeList,tBL,tST,numTapers,1,saveDir);
else   pwrST= nan(5,5,5,5,5); pwrBL= nan(5,5,5,5,5);  % dummy
end

% 3. some updates to parameter combinations
parameterCombinations.sValsUnique(parameterCombinations.sValsUnique<0 | parameterCombinations.sValsUnique>11) = 11.5; % for FullScreen, use same val
[~,indor] = intersect(parameterCombinations.oValsUnique,[22 67 112 157]); % for oris at gap of 22.5 
parameterCombinations.oValsUnique(indor) = parameterCombinations.oValsUnique(indor)+0.5;

end

function [powerST,powerBL,freqVals,parameterCombinations] = getPowerMT(folderSourceString,gridType,subjectName,expDate,protocolName,electrodeList,timeRangeBL,timeRangeST,numTapers,saveFlag,folderSave)

numElectrodes = length(electrodeList);
if numTapers>1,  tapers = [ceil((numTapers+1)/2) numTapers]; 
else             tapers = [1 1];
end 

fileSave = fullfile(folderSave,[subjectName expDate protocolName gridType '_N' num2str(numElectrodes) '_st_' num2str(1000*timeRangeST(1)) '_' ...
    num2str(1000*timeRangeST(2)) '_bl_' num2str(1000*timeRangeBL(1)) '_' num2str(1000*timeRangeBL(2)) '_tapers_' num2str(tapers(1)) '_' num2str(tapers(2)) '.mat']);

if exist(fileSave,'file')
    disp(['Loading ' fileSave]);
    x=load(fileSave);
    if isequal(x.electrodeList,electrodeList)
        powerST = x.powerST;
        powerBL = x.powerBL;
        freqVals= x.freqVals;
        parameterCombinations= x.parameterCombinations;
    end
else
    disp(['Generating ' fileSave]);
    
    % Get bad trials
    load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','badTrials.mat'));
    % Get timeVals
    load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','LFP','lfpInfo.mat')); 

    % Multi-taper
    Fs  =  round(1/(timeVals(2)-timeVals(1)));
    params.tapers   = tapers;
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = [0 250];
    params.trialave = 0;

    stPos  = find(timeVals>=timeRangeST(1),1) + (0:round(Fs*diff(timeRangeST))-1);
    blPos  = find(timeVals>=timeRangeBL(1),1) + (0:round(Fs*diff(timeRangeBL))-1);

     %%%%%%%%%%%%%%%%%%%%%%%%%%% Get Good StimPos %%%%%%%%%%%%%%%%%%%%%%%%%%
    clear goodStimPos
    parameterCombinations = load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData','parameterCombinations.mat')); % ParameterCombinations
    
    aPos=1;ePos=1;          % azimuth, elevation
    numSizes = length(parameterCombinations.sValsUnique);
    numSFs   = length(parameterCombinations.fValsUnique);             % spatial freqs
    numCons  = length(parameterCombinations.cValsUnique);
    numOris  = length(parameterCombinations.oValsUnique);    % orientations
    numTFs   = length(parameterCombinations.tValsUnique);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:numElectrodes
        electrodeNum = electrodeList(i,1);
             % Get LFP Data
        load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','LFP',['elec' num2str(electrodeNum) '.mat']));

        for sPos = 1:numSizes
            for fPos=1:numSFs 
                for cPos = 1:numCons
                    for oPos = 1:numOris
                        for tPos = 1:numTFs
                            clear goodPos
                            goodPos = setdiff(parameterCombinations.parameterCombinations{aPos,ePos,sPos,fPos,oPos,cPos,tPos},badTrials);
                            [powerST{i,sPos,fPos,cPos,oPos,tPos}]         = mtspectrumc(analogData(goodPos,stPos)',params);
                            [powerBL{i,sPos,fPos,cPos,oPos,tPos},freqVals]= mtspectrumc(analogData(goodPos,blPos)',params);
                            
                        end
                    end 
                end   
            end 
        end 
    end
    if saveFlag
        save(fileSave,'powerST','powerBL','freqVals','electrodeList','parameterCombinations');
    end  
end 
end 

function [pwrChangeInBand,meanPsdST,meanPsdBL,bandPwrST,bandPwrBL] = getPowerChanges(powerST,powerBL,freqVals,fBands,fBad)
% function to get change in power in given ranges. fBands is {[range1] [range2]}

[numElecs,numStimuli1,numStimuli2,numStimuli3,numStimuli4,numStimuli5] = size(powerBL);

meanPsdBL    = zeros([numElecs,numStimuli1,numStimuli2,numStimuli3,numStimuli4,numStimuli5,length(freqVals)]);
meanPsdST    = zeros([numElecs,numStimuli1,numStimuli2,numStimuli3,numStimuli4,numStimuli5,length(freqVals)]);
meanPsdBLall = nan(numElecs,length(freqVals));
bandPwrST    = zeros(numElecs,numStimuli1,numStimuli2,numStimuli3,numStimuli4,numStimuli5,length(fBands));
bandPwrBL    = zeros(numElecs,numStimuli1,numStimuli2,numStimuli3,numStimuli4,numStimuli5,length(fBands));
pwrChangeInBand = zeros(numElecs,numStimuli1,numStimuli2,numStimuli3,numStimuli4,numStimuli5,length(fBands));
[~,fBadInds] = intersect(freqVals,fBad);
    
for elec = 1:numElecs
    for pos1 = 1:numStimuli1  
        for pos2 = 1:numStimuli2       
            for pos3 = 1:numStimuli3   
                for pos4 = 1:numStimuli4   
                    for pos5 = 1:numStimuli5 
                        allLFPDataBL = (powerBL{elec,pos1,pos2,pos3,pos4,pos5});
                        allLFPDataST = (powerST{elec,pos1,pos2,pos3,pos4,pos5});
                        meanPsdBL(elec,pos1,pos2,pos3,pos4,pos5,:) = mean(allLFPDataBL,2);          % mean across trials
                        meanPsdST(elec,pos1,pos2,pos3,pos4,pos5,:) = mean(allLFPDataST,2);  
                        for ff = 1:length(fBands)
                            f_inds = freqVals>=fBands{ff}(1) & freqVals<=fBands{ff}(2);
                            f_inds(fBadInds) = false;
                            bandPwrST(elec,pos1,pos2,pos3,pos4,pos5,ff) = sum(squeeze(meanPsdST(elec,pos1,pos2,pos3,pos4,pos5,f_inds))); 
                            bandPwrBL(elec,pos1,pos2,pos3,pos4,pos5,ff) = sum(squeeze(meanPsdBL(elec,pos1,pos2,pos3,pos4,pos5,f_inds))); 
                        end 
                    end
                end
            end 
        end 
    end % end stimulus
end

meanPsdBLall(1:numElecs,:) = squeeze(nanmean(nanmean(nanmean(nanmean(nanmean(meanPsdBL,6),5),4),3),2));        % get average baseline across all stimuli
meanPsdBLall2 = permute(repmat(meanPsdBLall,[1 1 numStimuli1 numStimuli2 numStimuli3 numStimuli4 numStimuli5]),[1 3 4 5 6 7 2]); % get similar to meanPsdST
deltaPsd = meanPsdST./meanPsdBL;    

for ff = 1:length(fBands)
    f_inds = freqVals>=fBands{ff}(1) & freqVals<=fBands{ff}(2);
    f_inds(fBadInds)= false;
    pwrChangeInBand(:,:,:,:,:,:,ff) = (bandPwrST(:,:,:,:,:,:,ff) ./ bandPwrBL(:,:,:,:,:,:,ff));
end

end
