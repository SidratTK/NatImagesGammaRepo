% This function takes the hue patch predictions (allStimParams) and actual
% power values (allPower) for different images. It computes the predicted
% power when one or more of the patch properties are used for prediction.

% The predictions are for the following conditions:
% 1. Only H. Remaining to their max values 
% 2. Only S.
% 3. Only V.
% 4. Only R.
% 5. HS
% 6. HSV
% 7. HSVR (full model)

% The predicted power values for the full model are also returned.

% Correlation is computed for two cases:
% 1. All images (correlationValsFull)
% 2. selected images for which predition is deemed non-trivial (r>rCutoff)
% (correlationValsSelected)

function [correlationValsFull, correlationValsSelected, predictionString, predictedPower, selectedImageIndices] = getAllCorrelations(subjectName,allStimParams,allPower,rCutoff)

if ~exist('rCutoff','var');         rCutoff = 0;                      end
rMax = 10; % Large radius for which the gamma vs radius function saturates

numStimuli = length(allStimParams);
selectedImageIndices = [];
for i=1:numStimuli
    if (allStimParams{i}.radiusDeg > rCutoff) && ~(allStimParams{i}.gaborParams.categoryGabor)
        selectedImageIndices = cat(2,selectedImageIndices,i);
    end
end

% Only use Hue
tmpStimParams = allStimParams;
for i=1:numStimuli
    tmp = tmpStimParams{i};
    tmp.sat = 1;
    tmp.contrastPC = 100;
    tmp.spatialFreqPhaseDeg = 90;
    tmp.radiusDeg = rMax;
    tmpStimParams{i} = tmp;
end
predictionString{1} = 'H';
[correlationValsFull(1),correlationValsSelected(1)] = getCorrelations(subjectName,tmpStimParams,allPower,selectedImageIndices);

% Only use Sat
tmpStimParams = allStimParams;
for i=1:numStimuli
    tmp = tmpStimParams{i};
    tmp.hueDeg = 0; % Choose red
    tmp.contrastPC = 100;
    tmp.spatialFreqPhaseDeg = 90;
    tmp.radiusDeg = rMax;
    tmpStimParams{i} = tmp;
end
predictionString{2} = 'S';
[correlationValsFull(2),correlationValsSelected(2)] = getCorrelations(subjectName,tmpStimParams,allPower,selectedImageIndices);

% Only use Val - leave both contrastPC and spatialFreqPhaseDeg unchanged
tmpStimParams = allStimParams;
for i=1:numStimuli
    tmp = tmpStimParams{i};
    tmp.hueDeg = 0; % Choose red
    tmp.sat = 1;
    tmp.radiusDeg = rMax;
    tmpStimParams{i} = tmp;
end
predictionString{3} = 'V';
[correlationValsFull(3),correlationValsSelected(3)] = getCorrelations(subjectName,tmpStimParams,allPower,selectedImageIndices);

% Only use Size
tmpStimParams = allStimParams;
for i=1:numStimuli
    tmp = tmpStimParams{i};
    tmp.hueDeg = 0; % Choose red
    tmp.sat = 1;
    tmp.contrastPC = 100;
    tmp.spatialFreqPhaseDeg = 90;
    tmpStimParams{i} = tmp;
end
predictionString{4} = 'R';
[correlationValsFull(4),correlationValsSelected(4)] = getCorrelations(subjectName,tmpStimParams,allPower,selectedImageIndices);

% Only Hue and Sat
tmpStimParams = allStimParams;
for i=1:numStimuli
    tmp = tmpStimParams{i};
    tmp.contrastPC = 100;
    tmp.spatialFreqPhaseDeg = 90;
    tmp.radiusDeg = rMax;
    tmpStimParams{i} = tmp;
end
predictionString{5} = 'HS';
[correlationValsFull(5),correlationValsSelected(5)] = getCorrelations(subjectName,tmpStimParams,allPower,selectedImageIndices);

% Use HSV
tmpStimParams = allStimParams;
for i=1:numStimuli
    tmp = tmpStimParams{i};
    tmp.radiusDeg = rMax;
    tmpStimParams{i} = tmp;
end
predictionString{6} = 'HSV';
[correlationValsFull(6),correlationValsSelected(6)] = getCorrelations(subjectName,tmpStimParams,allPower,selectedImageIndices);

% Use HSVR (full model)
predictionString{7} = 'HSVR';
[correlationValsFull(7),correlationValsSelected(7),predictedPower] = getCorrelations(subjectName,allStimParams,allPower,selectedImageIndices); % Full model
    
end

function [rFull,rSelected,predictedPower] = getCorrelations(subjectName,tmpStimParams,allPower,selectedImageIndices)
numStimuli = length(tmpStimParams);

predictedPower = zeros(1,numStimuli);
for i=1:numStimuli
    predictedPower(i) = getPredictedGamma(subjectName,tmpStimParams{i});
end
tmp = corrcoef(allPower,predictedPower);
rFull = tmp(1,2);
if length(selectedImageIndices)>2 % Need to have at least 3 data points. Otherwise correlations are trivially 1 or -1.
    tmp = corrcoef(allPower(selectedImageIndices),predictedPower(selectedImageIndices));
    rSelected = tmp(1,2);
else
    rSelected = 0;
end
end