% get Gamma Parameters
% reads data files and fits curves to data to model gamma.
% inputs: subjectName = which monkey? 'alpaH','kesariH'
%         tBL,tST = baseline & stimlus time periods
%         fBandGamma = frequency range of interest
%         fBad = drop these frequencies from power calculation
%         HueFlag = read hue protocols(1) or grating protocols(0)
% outputs:modelParams = median across electrodes of parameters
%         modelParamsAll = values for all electrodes
%         modelParamsFold= all electrodes and their median parameters in 2 folds
%         gammaVals = struct containg values of power responses & electrodes protocol wise

function [modelParams,modelParamsAll,modelParamsFold,gammaVals]= getGammaParameters(subjectName,tBL,tST,fBandGamma,fBad,HueFlag)
gridType='Microelectrode';
folderSourceString=pwd; % where 'data' folder with extracted data is kept

if HueFlag % hue parameters
    [modelParams,modelParamsAll,modelParamsFold,gammaVals]= getGammaParametersHue(subjectName,tBL,tST,fBandGamma,fBad);
else       % grating parameters
    [modelParams,modelParamsAll,modelParamsFold,gammaVals]= getGammaParametersGrating(subjectName,tBL,tST,fBandGamma,fBad,gridType,folderSourceString);
end
end

function [modelParams,modelParamsAll,modelParamsFold,gammaVals]= getGammaParametersGrating(subjectName,tBL,tST,fBandGamma,fBad,gridType,folderSourceString)

[modelParams.sfCenter, modelParams.sfSigma,          ~,  modelParamsAllSf,  coeffNamesSf,  modelParamsFoldSf]  = getSfParams(subjectName,gridType,folderSourceString,tBL,tST,fBandGamma,fBad);
[modelParams.oriCenter,modelParams.oriSpread,    sfOri,  modelParamsAllOri, coeffNamesOri, modelParamsFoldOri] = getOrientationParams(subjectName,gridType,folderSourceString,tBL,tST,fBandGamma,fBad);
[modelParams.sizeSlope,modelParams.sizeAtHalfMax,sizeOri,modelParamsAllSize,coeffNamesSize,modelParamsFoldSize]= getSizeParams(subjectName,gridType,folderSourceString,tBL,tST,fBandGamma,fBad);
[modelParams.conSlope, modelParams.conAtHalfMax, conOri, modelParamsAllCon, coeffNamesCon, modelParamsFoldCon] = getContrastParams(subjectName,gridType,folderSourceString,tBL,tST,fBandGamma,fBad);

gammaVals.sfOri    = sfOri;
gammaVals.sizeOri  = sizeOri;
gammaVals.conOri   = conOri;
modelParamsAll(1:2)= modelParamsAllSf; modelParamsAll(3:4)=modelParamsAllOri; modelParamsAll(5:6)=modelParamsAllSize; modelParamsAll(7:8)=modelParamsAllCon;
modelParamsAll{9}  = [coeffNamesSf, coeffNamesOri, coeffNamesSize, coeffNamesCon];% [coeffNamesSfOri,coeffNamesSize,coeffNamesCon];%  
modelParamsFold.all= cat(2,modelParamsFoldSf.all, modelParamsFoldOri.all, modelParamsFoldSize.all, modelParamsFoldCon.all);
modelParamsFold.median=cat(2,modelParamsFoldSf.median, modelParamsFoldOri.median, modelParamsFoldSize.median, modelParamsFoldCon.median);

disp('calculated grating params');
end

function [modelParams,modelParamsAll,modelParamsFold,gammaVals]= getGammaParametersHue(subjectName,tBL,tST,fBandGamma,fBad)
if strcmp(subjectName(end),'H')          % if alpaH, kesariH
    subjectName = subjectName(1:end-1);   % use alpa, kesari
end
[modelParams.centerR,modelParams.spreadR,modelParams.centerC,modelParams.spreadC,modelParams.gainCbyR,...
    hueScreen,modelParamsAllHue,coeffNamesHue,modelParamsFoldHue]= getHueParams(subjectName,tBL,tST,fBandGamma,fBad);
[modelParams.sizeSlope,modelParams.sizeAtHalfMax,hueSize,modelParamsAllSize,coeffNamesSize,modelParamsFoldSize]= getHueSizeParams(subjectName,tBL,tST,fBandGamma,fBad);
[modelParams.satSlope,modelParams.satAtHalfMax,hueSatFullscreen,modelParamsAllSat,coeffNamesSat,modelParamsFoldSat]  = getHueSaturationParams(subjectName,tBL,tST,fBandGamma,fBad);
[modelParams.valSlope,modelParams.valIntercept,hueValFullscreen,modelParamsAllVal,coeffNamesVal,modelParamsFoldVal]  = getHueValueParams(subjectName,tBL,tST,fBandGamma,fBad);

gammaVals.hueScreen        = hueScreen;
gammaVals.hueSize          = hueSize;
gammaVals.hueSatFullscreen = hueSatFullscreen;
gammaVals.hueValFullscreen = hueValFullscreen;
modelParamsAll(1:5) = modelParamsAllHue; modelParamsAll(6:7)=modelParamsAllSize; modelParamsAll(8:9)=modelParamsAllSat; modelParamsAll(10:11)=modelParamsAllVal;
modelParamsAll{12}  = [coeffNamesHue, coeffNamesSize, coeffNamesSat, coeffNamesVal];
modelParamsFold.all = cat(2,modelParamsFoldHue.all, modelParamsFoldSize.all, modelParamsFoldSat.all, modelParamsFoldVal.all);
modelParamsFold.median=cat(2,modelParamsFoldHue.median, modelParamsFoldSize.median, modelParamsFoldSat.median, modelParamsFoldVal.median);

disp('calculated hue params');
end

%%%%
function [sfCenterMedian,sfSigmaMedian,sfOriVals,modelParamsAll,coeffNames,modelParamsFold] = getSfParams(subjectName,gridType,folderSourceString,tBL,tST,fBandGamma,fBad)
% define sf function
[~,gaussianFunctionStr] = gaussianFunctionOfLogx();
gaussianFunctionWithGainStr = [gaussianFunctionStr,'.*coeff(3) + coeff(4)'];
gaussianFunctionWithGain    = str2func(gaussianFunctionWithGainStr);
% get sf data 
datadir = 'savedData';
[expDates, protocolNames, ~, electrodeList, ~] = getProtocolInfoGratings(subjectName,gridType,folderSourceString,'SFOri',datadir);
[delGamma,~,~,~,~,parameterCombinations,~,delGammaF]= getEnergy(subjectName,expDates{1}{1},protocolNames{1}{1},fBandGamma,fBad,tBL,tST,electrodeList{1});
delGamma = squeeze(delGamma);
delGammaF= cellfun(@(x) squeeze(x), delGammaF,'un',0);
useSFs   = repmat(parameterCombinations.fValsUnique',[1 size(delGamma,3)]); % parameterCombinations.fValsUnique;%
% get estimates for params
for el = 1:length(electrodeList{1})
    delGammaUse = squeeze(delGamma(el,:,:)); %squeeze(mean(delGamma(el,:,:),3)); % 
    [coeffs,coeffNames]  = optimiseGaussianCoeffs(delGammaUse(:), useSFs(:), gaussianFunctionWithGain);
    modelParamsAll{1}(el,1) = coeffs(1); modelParamsAll{2}(el,1)=coeffs(2);
    for ff=1:2 % 2 folds
        delGammaUse = squeeze(delGammaF{ff}(el,:,:));% squeeze(mean(delGammaF{ff}(el,:,:),3)); 
        [coeffs,coeffNames]  = optimiseGaussianCoeffs(delGammaUse(:), useSFs(:), gaussianFunctionWithGain);
        modelParamsFold.all{ff,1}(el,1) = coeffs(1); modelParamsFold.all{ff,2}(el,1)=coeffs(2); 
    end
end
coeffNames = coeffNames(1:2);
% get median across elecs
sfCenterMedian=median(modelParamsAll{1}); sfSigmaMedian=median(modelParamsAll{2});
modelParamsFold.median=cellfun(@(x) median(x),modelParamsFold.all,'un',0);
sfOriVals.delGammaFold = delGammaF;
sfOriVals.delGamma = delGamma;
sfOriVals.spatialFreqCPD = parameterCombinations.fValsUnique; 
sfOriVals.orientationDeg = parameterCombinations.oValsUnique;
sfOriVals.radiusDeg  = 3*parameterCombinations.sValsUnique;  % radius/sigma= 3 in display s/w KNOT. doesnt matter for Full Screen
sfOriVals.contrastPC = parameterCombinations.cValsUnique; 
sfOriVals.electrodes = electrodeList{1};
end
function [oriCenterMedian,oriSpreadMedian,sfOriVals,modelParamsAll,coeffNames,modelParamsFold]= getOrientationParams(subjectName,gridType,folderSourceString,tBL,tST,fBandGamma,fBad)
% define ori function
[~,vmFunctionStr] = vmFunctionOf2Radx();
vmFunctionWithGainStr = [vmFunctionStr,'.*coeff(3) + coeff(4)'];
vmFunctionWithGain    = str2func(vmFunctionWithGainStr);
% get ori data 
datadir = 'savedData'; 
[expDates, protocolNames, ~, electrodeList, ~] = getProtocolInfoGratings(subjectName,gridType,folderSourceString,'SFOri',datadir);
[delGamma,~,~,~,~,parameterCombinations,~,delGammaF]= getEnergy(subjectName,expDates{1}{1},protocolNames{1}{1},fBandGamma,fBad,tBL,tST,electrodeList{1});
delGamma = squeeze(delGamma);
delGammaF= cellfun(@(x) squeeze(x), delGammaF,'un',0);
useOris  = repmat(parameterCombinations.oValsUnique,[size(delGamma,2) 1]); % parameterCombinations.oValsUnique;%
% get estimates for params
for el = 1:length(electrodeList{1})
    delGammaUse = squeeze(delGamma(el,:,:)); %squeeze(mean(delGamma(el,:,:),2)); %
    [coeffs,coeffNames]  = optimiseVmCoeffs(delGammaUse(:), useOris(:), vmFunctionWithGain);
    modelParamsAll{1}(el,1) = coeffs(1); modelParamsAll{2}(el,1)=coeffs(2);
    for ff=1:2
        delGammaUse = squeeze(delGammaF{ff}(el,:,:)); %squeeze(mean(delGammaF{ff}(el,:,:),2)); 
        [coeffs,coeffNames]  = optimiseVmCoeffs(delGammaUse(:), useOris(:), vmFunctionWithGain);
        modelParamsFold.all{ff,1}(el,1) = coeffs(1); modelParamsFold.all{ff,2}(el,1) = coeffs(2); 
    end
end
coeffNames = coeffNames(1:2);
% get median across elecs
oriCenterMedian = rad2deg(circ_median(deg2rad(modelParamsAll{1}))); oriSpreadMedian = median(modelParamsAll{2});
oriCenterMedian(oriCenterMedian<0) = oriCenterMedian(oriCenterMedian<0)+360; 
for ff=1:2
    modelParamsFold.median{ff,1} = rad2deg(circ_median(deg2rad(modelParamsFold.all{ff,1})));
    modelParamsFold.median{ff,1}(modelParamsFold.median{ff,1}<0)=modelParamsFold.median{ff,1}(modelParamsFold.median{ff,1}<0)+360;
    modelParamsFold.median{ff,2} = median(modelParamsFold.all{ff,2});
end 
sfOriVals.delGammaFold = delGammaF;
sfOriVals.delGamma = delGamma;
sfOriVals.spatialFreqCPD = parameterCombinations.fValsUnique; 
sfOriVals.orientationDeg = parameterCombinations.oValsUnique;
sfOriVals.radiusDeg  = 3*parameterCombinations.sValsUnique;  % radius/sigma= 3 in display s/w KNOT. doesnt matter for Full Screen
sfOriVals.contrastPC = parameterCombinations.cValsUnique; 
sfOriVals.electrodes = electrodeList{1};
end
function [sizeSlopeMedian,sizeAtHalfMaxMedian,sizeOriVals,modelParamsAll,coeffNames,modelParamsFold] = getSizeParams(subjectName,gridType,folderSourceString,tBL,tST,fBandGamma,fBad)
% define size function
[~,sigmoidFunctionStr]           = sigmoidalFunctionOfLogx();
sigmoidFunctionofLogxWithGainStr = [sigmoidFunctionStr,'.*coeff(3) + coeff(4)'];
sigmoidFunctionofLogxWithGain    = str2func(sigmoidFunctionofLogxWithGainStr);
% get size data electrode wise
datadir = 'savedData';
[expDates, protocolNames, ~, electrodeList, ~] = getProtocolInfoGratings(subjectName,gridType,folderSourceString,'SizeOri',datadir);
for el = 1:length(electrodeList{1})
    [delPwrBand,~,~,~,~,parameterCombinations,~,delPwrBandF]= getEnergy(subjectName,expDates{1}{el},protocolNames{1}{el},fBandGamma,fBad,tBL,tST,electrodeList{1}(el));
    sInd = parameterCombinations.sValsUnique>0.09;   % some alpaH protocols have an extra size. take the ones above 0.09
    parameterCombinations.sValsUnique = parameterCombinations.sValsUnique(sInd);  % are same for different days - checked
    delGamma(el,:,:) = squeeze(delPwrBand(:,sInd,:,:,:,:));
    delGammaUse = squeeze(delPwrBand(:,sInd,:,:,:,:));%mean(squeeze(delPwrBand(:,sInd,:,:,:,:)),2);  %
    useSizes= repmat(3*parameterCombinations.sValsUnique',[1 size(delGammaUse,2)]); %3*parameterCombinations.sValsUnique;% 3* because rad/sig=3 for size protocols
    % get estimates for params
    [coeffs,coeffNames]  = optimiseSigmoidalCoeffs(delGammaUse(:), useSizes(:), sigmoidFunctionofLogxWithGain);
    modelParamsAll{1}(el,1) = coeffs(1); modelParamsAll{2}(el,1)=coeffs(2); 
    for ff = 1:2 % folds
        delGammaF{ff}(el,:,:) = squeeze(delPwrBandF{ff}(:,sInd,:,:,:));
        delGammaUse = squeeze(delPwrBandF{ff}(:,sInd,:,:,:)); %  mean(squeeze(delPwrBandF{ff}(:,sInd,:,:,:)),2); %
        [coeffs,~]  = optimiseSigmoidalCoeffs(delGammaUse(:), useSizes(:), sigmoidFunctionofLogxWithGain);
        modelParamsFold.all{ff,1}(el,1) = coeffs(1); modelParamsFold.all{ff,2}(el,1)=coeffs(2); 
    end
end
coeffNames = cellfun(@(x) cat(2,'size ',x),coeffNames(1:2),'un',0);
% get median across elecs
sizeSlopeMedian=median(modelParamsAll{1}); sizeAtHalfMaxMedian=median(modelParamsAll{2});
modelParamsFold.median=cellfun(@(x) median(x),modelParamsFold.all,'un',0);
sizeOriVals.delGammaFold = delGammaF;
sizeOriVals.delGamma = delGamma;
sizeOriVals.spatialFreqCPD = parameterCombinations.fValsUnique; 
sizeOriVals.orientationDeg = parameterCombinations.oValsUnique; 
sizeOriVals.radiusDeg  = 3*parameterCombinations.sValsUnique;  % radius by sigma is 3 in display software KNOT
sizeOriVals.contrastPC = parameterCombinations.cValsUnique; 
sizeOriVals.electrodes = electrodeList{1}; 
end
function [conSlopeMedian,conAtHalfMaxMedian,conOriVals,modelParamsAll,coeffNames,modelParamsFold] = getContrastParams(subjectName,gridType,folderSourceString,tBL,tST,fBandGamma,fBad)
% define con function
[~,sigmoidFunctionStr]     = sigmoidalFunctionOfLogx();
sigmoidFunctionOfLogxWithGainStr = [sigmoidFunctionStr,'.*coeff(3) + coeff(4)'];
sigmoidFunctionOfLogxWithGain    = str2func(sigmoidFunctionOfLogxWithGainStr);
% get con data 
datadir = 'savedData'; % for data info
[expDates, protocolNames, ~, electrodeList, ~] = getProtocolInfoGratings(subjectName,gridType,folderSourceString,'ConOri',datadir);
[delPwrBand,~,~,~,~,parameterCombinations,~,delPwrBandF]= getEnergy(subjectName,expDates{1}{1},protocolNames{1}{1},fBandGamma,fBad,tBL,tST,electrodeList{1});

% ensure these are full screen protocols and not flickering
sInd = parameterCombinations.sValsUnique==11.5; tInd = parameterCombinations.tValsUnique==0;
parameterCombinations.sValsUnique = parameterCombinations.sValsUnique(sInd);
parameterCombinations.tValsUnique = parameterCombinations.tValsUnique(tInd);
delGamma = squeeze(delPwrBand(:,sInd,:,:,:,tInd));    %  as conori
delGammaF= cellfun(@(x) squeeze(x(:,sInd,:,:,:,tInd)),delPwrBandF,'un',0);
useCns   = repmat(parameterCombinations.cValsUnique',[1 size(delGamma,3)]);% parameterCombinations.cValsUnique;%   
% get estimates for params
for el = 1:length(electrodeList{1})
    delGammaUse = squeeze(delGamma(el,:,:)); %squeeze(mean(delGamma(el,:,:),3)); %
    [coeffs,coeffNames]  = optimiseSigmoidalCoeffs(delGammaUse(:), useCns(:), sigmoidFunctionOfLogxWithGain);
    modelParamsAll{1}(el,1) = coeffs(1); modelParamsAll{2}(el,1)=coeffs(2);
    for ff=1:2
        delGammaUse = squeeze(delGammaF{ff}(el,:,:)); %squeeze(mean(delGammaF{ff}(el,:,:),3)); %
        [coeffs,coeffNames]  = optimiseSigmoidalCoeffs(delGammaUse(:), useCns(:), sigmoidFunctionOfLogxWithGain);
        modelParamsFold.all{ff,1}(el,1) = coeffs(1); modelParamsFold.all{ff,2}(el,1) = coeffs(2); 
    end
end
coeffNames = cellfun(@(x) cat(2,'con ',x),coeffNames(1:2),'un',0);
% get median across elecs
conSlopeMedian=median(modelParamsAll{1}); conAtHalfMaxMedian=median(modelParamsAll{2});
modelParamsFold.median=cellfun(@(x) median(x),modelParamsFold.all,'un',0);
conOriVals.delGammaFold = delGammaF;
conOriVals.delGamma = delGamma;
conOriVals.spatialFreqCPD = parameterCombinations.fValsUnique; 
conOriVals.orientationDeg = parameterCombinations.oValsUnique;
conOriVals.radiusDeg  = 3*parameterCombinations.sValsUnique;  % radius by sigma is 3 in display software KNOT
conOriVals.contrastPC = parameterCombinations.cValsUnique; 
conOriVals.electrodes = electrodeList{1};
end
function [sizeSlopeMedian,sizeAtHalfMaxMedian,hueSizeVals,modelParamsAll,coeffNames,modelParamsFold] = getHueSizeParams(subjectName,tBL,tST,fBandGamma,fBad)
% define size function
[~,sigmoidFunctionStr]     = sigmoidalFunctionOfLogx();
sigmoidFunctionofLogxWithGainStr = [sigmoidFunctionStr,'.*coeff(3) + coeff(4)'];
sigmoidFunctionofLogxWithGain    = str2func(sigmoidFunctionofLogxWithGainStr);
% get size data electrode wise - all elecs on one day
datadir = 'savedData';
[expDates, protocolNames, ~, electrodeList, ~] = getProtocolInfoHues(subjectName,'HueSize',datadir);
if strcmp(subjectName,'kesari'),rbys=1; else, rbys=3; end  % this is input in knot
for el = 1:length(electrodeList{1})
    [delPwrBand,~,~,~,~,parameterCombinations,~,delPwrBandF]= getEnergy(subjectName,expDates{1}{1},protocolNames{1}{1},fBandGamma,fBad,tBL,tST,electrodeList{1}(el));
    useSizes1 = rbys*(parameterCombinations.sValsUnique); 
    delGammaUse = squeeze(delPwrBand(:,:,:,:,parameterCombinations.oValsUnique==0));
    delGamma(el,:) = delGammaUse; %   
    % get estimates for params
    [coeffs,coeffNames]  = optimiseSigmoidalCoeffs(delGammaUse(:), useSizes1(:), sigmoidFunctionofLogxWithGain);
    modelParamsAll{1}(el,1) = coeffs(1); modelParamsAll{2}(el,1)=coeffs(2); 
    for ff = 1:2 % folds
        delGammaUse = squeeze(delPwrBandF{ff}(:,:,:,:,parameterCombinations.oValsUnique==0));
        delGammaF{ff}(el,:) = delGammaUse; 
        [coeffs,~]  = optimiseSigmoidalCoeffs(delGammaUse(:), useSizes1(:), sigmoidFunctionofLogxWithGain);
        modelParamsFold.all{ff,1}(el,1) = coeffs(1); modelParamsFold.all{ff,2}(el,1)=coeffs(2); 
    end
end
coeffNames = cellfun(@(x) cat(2,'size ',x),coeffNames(1:2),'un',0);
% get median across elecs
sizeSlopeMedian=median(modelParamsAll{1}); sizeAtHalfMaxMedian=median(modelParamsAll{2});
modelParamsFold.median=cellfun(@(x) median(x),modelParamsFold.all,'un',0);
hueSizeVals.delGammaFold = delGammaF;
hueSizeVals.delGamma = delGamma;
hueSizeVals.sat    = parameterCombinations.fValsUnique; 
hueSizeVals.hueDeg = parameterCombinations.oValsUnique(parameterCombinations.oValsUnique==0); % only red used
hueSizeVals.radiusDeg = parameterCombinations.sValsUnique * rbys; 
hueSizeVals.val    = parameterCombinations.cValsUnique/100; % get 0-1 
hueSizeVals.electrodes= electrodeList{1}; 
end
function [rCenMedian,rSprMedian,cCenMedian,cSprMedian,gainCbyRmedian,hueScreenVals,modelParamsAll,coeffNames,modelParamsFold] = getHueParams(subjectName,tBL,tST,fBandGamma,fBad)
% define hue function
[~,vmFunctionStr]     = sum2VmFunctionsOfRadx();                
vmFunctionWithGainStr = [vmFunctionStr,'.*coeff(6) + coeff(7)']; 
vmFunctionWithGain    = str2func(vmFunctionWithGainStr);
datadir = 'savedData';
[expDates, protocolNames, protocolTypesAll, electrodeList, ~] = getProtocolInfoHues(subjectName,[],datadir);
% get  data of full screen protocols
[delPwrBand1,~,~,~,~,parameterCombinations,~,delPwrBand1F] = getEnergy(subjectName,expDates{strcmp(protocolTypesAll,'HueScreen')}{1},protocolNames{strcmp(protocolTypesAll,'HueScreen')}{1},fBandGamma,fBad,tBL,tST,electrodeList{strcmp(protocolTypesAll,'HueScreen')});
[delPwrBand2,~,~,~,~,parameterCombinations2,~,delPwrBand2F]= getEnergy(subjectName,expDates{strcmp(protocolTypesAll,'HueSatFullscreen')}{1},protocolNames{strcmp(protocolTypesAll,'HueSatFullscreen')}{1},fBandGamma,fBad,tBL,tST,electrodeList{strcmp(protocolTypesAll,'HueSatFullscreen')});
[delPwrBand3,~,~,~,~,parameterCombinations3,~,delPwrBand3F]= getEnergy(subjectName,expDates{strcmp(protocolTypesAll,'HueValFullscreen')}{1},protocolNames{strcmp(protocolTypesAll,'HueValFullscreen')}{1},fBandGamma,fBad,tBL,tST,electrodeList{strcmp(protocolTypesAll,'HueValFullscreen')});
% use huescreen & append data from hue-sat and hue-val protocols as well
useHues1 = parameterCombinations.oValsUnique;
useHues2 = parameterCombinations2.oValsUnique; 
useHues3 = parameterCombinations3.oValsUnique;
useHues  = (cat(1,useHues1(:),useHues2(:),useHues3(:))); %
delGamma1 = squeeze(delPwrBand1(:,1,1,1,:));   delGamma1F = cellfun(@(x) squeeze(x(:,1,1,1,:)),delPwrBand1F,'un',0); 
delGamma2 = squeeze(delPwrBand2(:,1,end,1,:)); delGamma2F = cellfun(@(x) squeeze(x(:,1,end,1,:)),delPwrBand2F,'un',0); % at max sat 
delGamma3 = squeeze(delPwrBand3(:,1,1,end,:)); delGamma3F = cellfun(@(x) squeeze(x(:,1,1,end,:)),delPwrBand3F,'un',0); % at max val 
for el = 1:length(electrodeList{1})
    delGammaUse = squeeze(cat(2,delGamma1(el,:)./delGamma1(el,1), delGamma2(el,:)./delGamma2(el,1), delGamma3(el,:)./delGamma3(el,1))); % divide by red response for each protocol.
    % get estimates for params
    [coeffs,coeffNames]  = optimiseSum2VmCoeffs(delGammaUse(:), useHues(:), vmFunctionWithGain);
    modelParamsAll{1}(el,1) = coeffs(1); modelParamsAll{2}(el,1) =coeffs(2); 
    modelParamsAll{3}(el,1) = coeffs(3); modelParamsAll{4}(el,1) =coeffs(4); 
    modelParamsAll{5}(el,1) = coeffs(5); 
    for ff=1:2
      delGammaUse = squeeze(cat(2,delGamma1F{ff}(el,:)./delGamma1F{ff}(el,1), delGamma2F{ff}(el,:)./delGamma2F{ff}(el,1), delGamma3F{ff}(el,:)./delGamma3F{ff}(el,1)));
      [coeffs,coeffNames]  = optimiseSum2VmCoeffs(delGammaUse(:), useHues(:), vmFunctionWithGain);
      modelParamsFold.all{ff,1}(el,1) = coeffs(1); modelParamsFold.all{ff,2}(el,1) =coeffs(2); 
      modelParamsFold.all{ff,3}(el,1) = coeffs(3); modelParamsFold.all{ff,4}(el,1) =coeffs(4); 
      modelParamsFold.all{ff,5}(el,1) = coeffs(5); 
    end
end
coeffNames = coeffNames(1:5);  
% get median across elecs
rCenMedian=circ_median(modelParamsAll{1}); rSprMedian=median(modelParamsAll{2});
cCenMedian=circ_median(modelParamsAll{3}); cSprMedian=median(modelParamsAll{4});
gainCbyRmedian= median(modelParamsAll{5});
rCenMedian(rCenMedian<0)= rCenMedian(rCenMedian<0)+2*pi;
cCenMedian(cCenMedian<0)= cCenMedian(cCenMedian<0)+2*pi;
modelParamsFold.median=cellfun(@(x) median(x),modelParamsFold.all,'un',0);
for ii=[1 3]  % ori center angles
    for ff=1:2
        modelParamsFold.median{ff,ii} = circ_median(modelParamsFold.all{ff,ii});
        modelParamsFold.median{ff,ii}(modelParamsFold.median{ff,ii}<0)=modelParamsFold.median{ff,ii}(modelParamsFold.median{ff,ii}<0)+2*pi;
    end
end
hueScreenVals.delGamma = delGamma1;
hueScreenVals.delGammaFold = delGamma1F;
hueScreenVals.hueDeg  = useHues1; 
hueScreenVals.radiusDeg = 3*parameterCombinations.sValsUnique; % because r/s = 3 in knot
hueScreenVals.sat = parameterCombinations.fValsUnique; 
hueScreenVals.val = parameterCombinations.cValsUnique/100; % get to 0-1; 
hueScreenVals.electrodes= electrodeList{1};
end
function [satSlopeMedian,satAtHalfMaxMedian,hueSatFullscreenVals,modelParamsAll,coeffNames,modelParamsFold] = getHueSaturationParams(subjectName,tBL,tST,fBandGamma,fBad)
% define sat function
[~,sigmoidFunctionStr]     = sigmoidalFunctionOfx();
sigmoidFunctionWithGainStr = [sigmoidFunctionStr,'.*coeff(3) + coeff(4)'];
sigmoidFunctionWithGain    = str2func(sigmoidFunctionWithGainStr);
datadir = 'savedData';
[expDates, protocolNames, ~, electrodeList, ~] = getProtocolInfoHues(subjectName,'HueSatFullscreen',datadir);
% get  data 
[delPwrBand,~,~,~,~,parameterCombinations,~,delPwrBandF]= getEnergy(subjectName,expDates{1}{1},protocolNames{1}{1},fBandGamma,fBad,tBL,tST,electrodeList{1});
delGamma = squeeze(delPwrBand); 
delGammaF= cellfun(@(x) squeeze(x), delPwrBandF, 'un',0);    % as sat ori
useSats = repmat(parameterCombinations.fValsUnique',[1 size(delGamma,3)]);
for el = 1:length(electrodeList{1})
    delGammaUse = squeeze(delGamma(el,:,:)); 
    % get estimates for params
    [coeffs,coeffNames]  = optimiseSigmoidalCoeffs(delGammaUse(:), useSats(:), sigmoidFunctionWithGain);
    modelParamsAll{1}(el,1) = coeffs(1); modelParamsAll{2}(el,1)=coeffs(2); 
    for ff=1:2
        delGammaUse = squeeze(delGammaF{ff}(el,:,:));   
        [coeffs,coeffNames]  = optimiseSigmoidalCoeffs(delGammaUse(:), useSats(:), sigmoidFunctionWithGain);
        modelParamsFold.all{ff,1}(el,1) = coeffs(1); modelParamsFold.all{ff,2}(el,1)=coeffs(2); 
    end 
end
coeffNames = cellfun(@(x) cat(2,'sat ',x),coeffNames(1:2),'un',0);
% get median across elecs
satSlopeMedian=median(modelParamsAll{1}); satAtHalfMaxMedian=median(modelParamsAll{2});
modelParamsFold.median=cellfun(@(x) median(x),modelParamsFold.all,'un',0);
hueSatFullscreenVals.delGammaFold = delGammaF;
hueSatFullscreenVals.delGamma = delGamma;
hueSatFullscreenVals.sat    = parameterCombinations.fValsUnique; 
hueSatFullscreenVals.hueDeg = parameterCombinations.oValsUnique;
hueSatFullscreenVals.radiusDeg = 3*parameterCombinations.sValsUnique; % because r/s = 3 in knot
hueSatFullscreenVals.val    = parameterCombinations.cValsUnique/100; % get to 0-1 
hueSatFullscreenVals.electrodes= electrodeList{1};
end
function [valSlopeMedian,valAtHalfMaxMedian,hueValFullscreenVals,modelParamsAll,coeffNames,modelParamsFold] = getHueValueParams(subjectName,tBL,tST,fBandGamma,fBad)
% define val function
[~,linearFunctionStr] = linearFunctionOfx(); 
linearFunctionx       = str2func(linearFunctionStr);   
datadir = 'savedData';
[expDates, protocolNames, ~, electrodeList, ~] = getProtocolInfoHues(subjectName,'HueValFullscreen',datadir);
% get  data 
[delPwrBand,~,~,~,~,parameterCombinations,~,delPwrBandF]= getEnergy(subjectName,expDates{1}{1},protocolNames{1}{1},fBandGamma,fBad,tBL,tST,electrodeList{1});
delGamma = squeeze(delPwrBand(:,:,:,parameterCombinations.cValsUnique>0,:)); % skip the first V val v=0
delGammaF= cellfun(@(x) squeeze(x(:,:,:,parameterCombinations.cValsUnique>0,:)), delPwrBandF, 'un',0);    % as val ori
useVals = parameterCombinations.cValsUnique(parameterCombinations.cValsUnique>0);  % drop the 0% contrst
useVals = repmat(useVals'/100,[1 size(delGamma,3)]);  %  this is a percentage. scale to 0-1
for el = 1:length(electrodeList{1})
    delGammaUse = squeeze(delGamma(el,:,:));
    delGammaUse = delGammaUse./max(delGammaUse(:));      % this to make max as 1. parameters learnt on scaled version
    % get estimates for params
    [coeffs,coeffNames]  = optimiseLinearCoeffs(delGammaUse(:), useVals(:), linearFunctionx);
    modelParamsAll{1}(el,1) = coeffs(1); modelParamsAll{2}(el,1)=coeffs(2); 
    for ff=1:2
        delGammaUse = squeeze(delGammaF{ff}(el,:,:)); % 
        delGammaUse = delGammaUse./max(delGammaUse(:));      % this to make max as 1. parameters learnt on scaled version
        [coeffs,coeffNames]  = optimiseLinearCoeffs(delGammaUse(:), useVals(:), linearFunctionx);
        modelParamsFold.all{ff,1}(el,1) = coeffs(1); modelParamsFold.all{ff,2}(el,1)=coeffs(2); 
    end 
end
coeffNames = cellfun(@(x) cat(2,'val ',x),coeffNames(1:2),'un',0);
% get median across elecs
valSlopeMedian=median(modelParamsAll{1}); valAtHalfMaxMedian=median(modelParamsAll{2});
modelParamsFold.median=cellfun(@(x) median(x),modelParamsFold.all,'un',0);
hueValFullscreenVals.delGammaFold = delGammaF;
hueValFullscreenVals.delGamma = delGamma;
hueValFullscreenVals.sat    = parameterCombinations.fValsUnique; 
hueValFullscreenVals.hueDeg = parameterCombinations.oValsUnique;
hueValFullscreenVals.radiusDeg = 3*parameterCombinations.sValsUnique; % because r/s = 3 in knot
hueValFullscreenVals.val    = useVals(:,1);%parameterCombinations.cValsUnique; %
hueValFullscreenVals.electrodes= electrodeList{1};
end
%%%%
function [y,gaussianFunctionStr] = gaussianFunctionOfLogx(x,coeff)
if ~exist('x','var'), x = 0; end
if ~exist('coeff','var'), coeff = [0 1]; end  % mu = coeff(1); sigma = coeff(2);
gaussianFunction = @(coeff,x) normpdf(log2(x),log2(coeff(1)),coeff(2))./normpdf(log2(coeff(1)),log2(coeff(1)),coeff(2)) ;  % scaled so max is 1
y = gaussianFunction(coeff,x);
gaussianFunctionStr = func2str(gaussianFunction);
end
function [y,vmFunctionStr] = vmFunctionOf2Radx(x,coeff)
if ~exist('x','var'), x = 0; end
if ~exist('coeff','var'), coeff = [0 1]; end  % thetaPreferred = coeff(1). concentration parameter =coeff(2)
vmFunction= @(coeff,x) circ_vmpdf(2*deg2rad(x),2*deg2rad(coeff(1)),coeff(2))./circ_vmpdf(2*deg2rad(coeff(1)),2*deg2rad(coeff(1)),coeff(2));  %  scaled so max is 1
y = vmFunction(coeff,x);
vmFunctionStr = func2str(vmFunction);
end
function [y,sigmoidFunctionStr] = sigmoidalFunctionOfLogx(x,coeff)
if ~exist('x','var'), x = 1; end
if ~exist('coeff','var'), coeff = [1 0]; end  % slope = coeff(1); offset = coeff(2);
sigmoidFunction = @(coeff,x) (1./(1+10.^(coeff(1)*(log2(coeff(2))-log2(x)))));  % put logCoef2 so that x and coeff are on the same scale / range
y = sigmoidFunction(coeff,x);                % y = 1./(1+10.^(slope*(x-offset)));
sigmoidFunctionStr = func2str(sigmoidFunction);
end
function [y,sigmoidFunctionStr] = sigmoidalFunctionOfx(x,coeff)
if ~exist('x','var'), x = 1; end
if ~exist('coeff','var'), coeff = [1 0]; end  % slope = coeff(1); offset = coeff(2);
sigmoidFunction = @(coeff,x) (1./(1+10.^(coeff(1)*(coeff(2)-x))));
y = sigmoidFunction(coeff,x);                % y = 1./(1+10.^-(slope*(x-offset)));
sigmoidFunctionStr = func2str(sigmoidFunction);
end
function [y,sum2VmFunctionStr] = sum2VmFunctionsOfRadx(x,coeff)
if ~exist('x','var'), x = 0; end
if ~exist('coeff','var'), coeff = [0 1 0 1 0]; end
sum2VmFunction =  @(coeff,x)( circ_vmpdf(deg2rad(x),coeff(1),coeff(2))./circ_vmpdf(coeff(1),coeff(1),coeff(2)) + ...  % keep max at 1.
                    coeff(5).*circ_vmpdf(deg2rad(x),coeff(3),coeff(4))./circ_vmpdf(coeff(3),coeff(3),coeff(4)));
y = sum2VmFunction(coeff,x);
sum2VmFunctionStr = func2str(sum2VmFunction);
end
function [y,linearFunctionStr] = linearFunctionOfx(x,coeff)
if ~exist('x','var'), x = 1; end
if ~exist('coeff','var'), coeff = [1 0]; end
linearFunction = @(coeff,x)(coeff(1).*x + coeff(2));%./(coeff(1).*100 + coeff(2)); % max of x=100. divide this to keep max y 1
y = linearFunction(coeff,x);
linearFunctionStr = func2str(linearFunction);
end

%%%%
function [fitResult,fitResultName] = optimiseGaussianCoeffs(yData, xData, gaussianFunctionWithGain)
% initial seed params:
[maxis,ind] = max(yData);  minis = min(yData);
sseed = var(log2(xData),yData);
X0 =    [ xData(ind)   sseed  maxis-minis  minis];  %
Lower = [ min(xData)   eps       0         0]; 
Upper = [ max(xData)   Inf      Inf       Inf];
fitResultName = {'sfCenter','sfSigma','gain','offset'};
% Set up fittype and options.
tols    = 1e-16;
options = optimset('Display','off','MaxFunEvals',5e2,'MaxIter',5e2,'TolX',tols,'TolFun',tols,'TolCon',tols );
% perform the fitting
[fitResult,~,~] = lsqcurvefit(gaussianFunctionWithGain,X0,xData,yData,Lower,Upper,options);
end
function [fitResult,fitResultName] = optimiseVmCoeffs(yData, xData, vmFunctionWithGain)
% initial seed params:
[maxis,ind] = max(yData);  minis = min(yData);
[~,kseed] = circ_var(2*deg2rad(xData),yData);
kseed = 1/kseed;
X0 =    [ xData(ind) kseed maxis-minis  minis];  %
Lower = [ 0          0     0            0]; 
Upper = [ 180        200   Inf          Inf];
fitResultName = {'oriCenter','oriSpread','gain','offset'};
% Set up fittype and options.
tols    = 1e-16;
options = optimset('Display','off','MaxFunEvals',5e2,'MaxIter',5e2,'TolX',tols,'TolFun',tols,'TolCon',tols );
% perform the fitting
[fitResult,~,~] = lsqcurvefit(vmFunctionWithGain,X0,xData,yData,Lower,Upper,options);
end
function [fitResult,fitResultName] = optimiseSigmoidalCoeffs(yData, xData, sigmoidFunctionWithGain)
% initial seed params:
xhalf = xData(2);    %  just a initial value
X0    = [1   xhalf                max(yData(:)) 0]; %   seeding for slope & half max & gain & offset
Lower = [0   min([0, min(xData)]) 0             0]; % 
Upper = [Inf max(xData)           Inf           Inf];   
fitResultName = {'slope','halfMax','gain','offset'};
% Set up fittype and options.
tols = 1e-16;
options = optimset('Display','off','MaxFunEvals',5e2,'MaxIter',5e2,'TolX',tols,'TolFun',tols,'TolCon',tols );
% perform the fitting
[fitResult,~,~] = lsqcurvefit(sigmoidFunctionWithGain,X0,xData,yData,Lower,Upper,options);
end
function [fitResult,fitResultName] = optimiseSum2VmCoeffs(yData, xData, sum3VmFunctionWithGain)
% initial seed params:
cen = [0 180]; %  keep vm centers at red(0), cyan(180)
cenR= deg2rad(cen);
minis = min(yData);
G(1) = mean(yData(xData==cen(1))-minis); % gains
G(2) = mean(yData(xData==cen(2))-minis)./G(1);   % norm by red gain
s = 5;         % start with a spread . 
cenlo= deg2rad([cen(1)-80 cen(2)-100]); % 
cenup= deg2rad([cen(1)+80 cen(2)+100]); %

X0 = [cenR(1)  s   cenR(2)  s    G(2)  G(1) minis ];%
Low= [cenlo(1) 0   cenlo(2) 0    0     0    0     ]; %
Upp= [cenup(1) 50  cenup(2) 50   Inf   Inf  Inf ];   % 
fitResultName{1} ='centerR'; fitResultName{2} ='spreadR';
fitResultName{3} ='centerC'; fitResultName{4} ='spreadC';
fitResultName{5} ='gainCbyR';
fitResultName{6} ='gain';    fitResultName{7}='offset';
% Set up fittype and options.
tols    = 1e-16;
options = optimset('Display','off','MaxFunEvals',5e2,'MaxIter',5e2,'TolX',tols,'TolFun',tols,'TolCon',tols );
% perform the fitting
[fitResult,~,~] = lsqcurvefit(sum3VmFunctionWithGain,X0,xData,yData,Low,Upp,options);
end
function [fitResult,fitResultName] = optimiseLinearCoeffs(yData, xData, linearFunctionofx)
% initial seed params:
sl = (yData(end)-yData(1))./(xData(end)-xData(1));
X0    = [sl  min(yData(:))];%   seeding for slope & intercept
Lower = [0   0            ];% 
Upper = [Inf Inf          ];   
fitResultName = {'slope','intercept'};
% Set up fittype and options.
tols = 1e-16;
options = optimset('Display','off','MaxFunEvals',5e2,'MaxIter',5e2,'TolX',tols,'TolFun',tols,'TolCon',tols );
% perform the fitting
[fitResult,~,~] = lsqcurvefit(linearFunctionofx,X0,xData,yData,Lower,Upper,options);
end
%%%%
function [expDates, protocolNames, protocolTypes, electrodes, elecType] = getProtocolInfoGratings(monkeyName,gridType,folderSourceString,protocolType,datadir)
if ~exist('gridType','var'),           gridType = 'Microelectrode';                      end
if ~exist('folderSourceString','var'), folderSourceString = pwd; end
if ~exist('datadir','var'),            datadir =  'savedData';  end
if ~exist('protocolType','var'),       protocolType = [];                                end
    
protocolTypes = {'SFOri', 'SizeOri', 'ConOri'};
expDates      = cell(length(protocolTypes),1);
protocolNames = cell(length(protocolTypes),1);
electrodes    = cell(length(protocolTypes),1);
elecType      = cell(length(protocolTypes),1);

[~,~,LFPElectrodeList,EcogElectrodeList,~] = getRFdetails(monkeyName,datadir); % dataDir has 'RFDetails.mat' file
for pp = [1 3]    % full screen protocols
    electrodes{pp} =  cat(2,LFPElectrodeList{1}',EcogElectrodeList{1}');  
    elecType{pp}   =  cat(2,repmat({'LFP'},[1 length(LFPElectrodeList{1})]), repmat({'ECoG'},[1 length(EcogElectrodeList{1})]));
end
        
if strcmp(monkeyName,'alpaH')
    expDates{1}      = {'210817' };  % SFOri FullScreen r/s=3
    protocolNames{1} = {'GRF_002'};
    expDates{3}      = {'290817' }; % conTFOri r/s=3
    protocolNames{3} = {'GRF_001'};
    electrodes{2}    =  [82      85      86      88     89];      % ecogs for size protocol
    elecType{2}      = {'ECoG', 'ECoG', 'ECoG', 'ECoG', 'ECoG'};
    expDates{2}      = {'130917' ,'080917' ,'080917' ,'170817' ,'260817'};    % size ori electrode wise r/s=3
    protocolNames{2} = {'GRF_001','GRF_001','GRF_002','GRF_005','GRF_002'};
    
elseif strcmp(monkeyName,'kesariH')
    expDates{1}      = {'270218' };  % SFOri FullScreen r/s=3 
    protocolNames{1} = {'GRF_001'};
    expDates{3}      = {'240118' };  % conOri r/s=3 
    protocolNames{3} = {'GRF_002'};
    electrodes{2}    = [85,        86,       88,       89,  ];    % ecogs for size protocol.
    elecType{2}      = {'ECoG',   'ECoG',   'ECoG',   'ECoG' };
    expDates{2}      = {'271217', '030118', '060118', '070118'};   % size ori electrode wise r/s=3 
    protocolNames{2} = {'GRF_002','GRF_002','GRF_002','GRF_001'};
    
elseif strcmp(monkeyName,'alpa')
    protocolTypes = {'SFOri'}; % only this
    expDates{1}   = {'230216'};
    protocolNames{1} = {'GRF_001'};

elseif strcmp(monkeyName,'kesari')
    protocolTypes = {'SFOri'};
    expDates{1}   = {'250416'};
    protocolNames{1} = {'GRF_001'};   

end

if isempty(protocolType) || strcmp(protocolType,'SizeOri')
    load(['savedData/elecListSizeOri',monkeyName,'.mat']);
    electrodes{2}   = cat(2, elecListAll, electrodes{2});
    elecType{2}     = cat(2, repmat({'LFP'},[1 length(elecListAll)]), elecType{2});
    expDates{2}     = cat(2, expDateAll, expDates{2});
    protocolNames{2}= cat(2, protocolNameAll, protocolNames{2});
end

if ~isempty(protocolType)
   pInd         = strcmp(protocolTypes,protocolType);
   expDates     = expDates(pInd);
   protocolNames= protocolNames(pInd);
   electrodes   = electrodes(pInd);
   elecType     = elecType(pInd);
end
end

function [expDates,protocolNames,protocolTypes,electrodes,elecType] = getProtocolInfoHues(monkeyName,protocolType,datadir)
if ~exist('datadir','var'),            datadir =  'savedData';  end
if ~exist('protocolType','var'),       protocolType = [];       end

protocolTypes = {'HueScreen', 'HueSize', 'HueSatFullscreen', 'HueValFullscreen','SatVal'};
expDates      = cell(length(protocolTypes),1);
protocolNames = cell(length(protocolTypes),1);
electrodes    = cell(length(protocolTypes),1);
elecType      = cell(length(protocolTypes),1);

[rfStatsDeg,~,LFPElectrodeList,EcogElectrodeList,~] = getRFdetails(monkeyName,datadir); % dataDir has 'RFDetails.mat' file
for pp = [1 3 4]   % full screen protocols
    electrodes{pp} =  cat(2,LFPElectrodeList{1}',EcogElectrodeList{1}');  
    elecType{pp}   =  cat(2,repmat({'LFP'},[1 length(LFPElectrodeList{1})]), repmat({'ECoG'},[1 length(EcogElectrodeList{1})]));
end

if strcmp(monkeyName,'alpaH')
    expDates{1}      = {'010917'};   % hue screen protocol
    protocolNames{1} = {'GRF_005'};
    
elseif strcmp(monkeyName,'kesariH')
    expDates{1}      = {'050318'};  % hue screen protocol
    protocolNames{1} = {'GRF_004'};
    
elseif strcmp(monkeyName,'alpa')
    expDates{1}      = {'301215'};  %
    protocolNames{1} = {'GRF_001'};
    expDates{2}      = {'140915' };
    protocolNames{2} = {'GRF_004'};  %  r/s=3
    stimAzi          = -2.05;  stimEle = -3.8; % stimulus location for size protocol
    dCutoff          = 0.3;                    % distance from stim center
    [electrodes{2}]  = getElectInStim(rfStatsDeg{1},electrodes{1},stimAzi,stimEle,dCutoff);
    elecType{2}      = repmat({'LFP'},[1 length(electrodes{2})]);
    expDates{3}      = {'021015'};    % 
    protocolNames{3} = {'GRF_001'} ;
    expDates{4}      = {'011015'};    % 
    protocolNames{4} = {'GRF_002'} ;
    expDates{5}      = {'140915'};    % 
    protocolNames{5} = {'GRF_003'} ;
    stimAzi          = -2.05;  stimEle   = -3.8; % stimulus location for SatVal
    dCutoff          = 0.5;        % distance from stim center
    [electrodes{5}]  = getElectInStim(rfStatsDeg{1},electrodes{1},stimAzi,stimEle,dCutoff);
    elecType{5}      =  repmat({'LFP'},[1 length(electrodes{5})]);
    
elseif strcmp(monkeyName,'kesari')
    expDates{1}      = {'050516'};  % 
    protocolNames{1} = {'GRF_001'};
    expDates{2}      = {'210716'};
    protocolNames{2} = {'GRF_006'};       % r/s=1
    stimAzi          = -0.8;  stimEle= -1.3; % stimulus location fr size protocol
    dCutoff          = 0.3;                 % distance from stim center
    [electrodes{2}]  = getElectInStim(rfStatsDeg{1},electrodes{1},stimAzi,stimEle,dCutoff);
    elecType{2}      = repmat({'LFP'},[1 length(electrodes{2})]);
    expDates{3}      = {'021015'};    % 
    expDates{3}      = {'250516'};    %
    protocolNames{3} = {'GRF_001'} ;
    expDates{4}      = {'270516'};    
    protocolNames{4} = {'GRF_001'} ;
end

if ~isempty(protocolType)
   pInd         = strcmp(protocolTypes,protocolType);
   expDates     = expDates(pInd);
   protocolNames= protocolNames(pInd);
   electrodes   = electrodes(pInd);
   elecType     = elecType(pInd);
end

end
function [useElecs, useElecInds] = getElectInStim(rfStats,electrodeList,a,e,dCutoff)
useElecs = []; useElecInds=[];
for j=1:length(electrodeList)
    azi = rfStats(electrodeList(j)).meanAzi;
    ele = rfStats(electrodeList(j)).meanEle;
    d   = sqrt(sum((azi-a)^2+(ele-e)^2)); % calc the distance between RF center and given location
    if d < dCutoff
        useElecs =   cat(2,useElecs,electrodeList(j));
        useElecInds= cat(2,useElecInds, j);
    end
end          
end

