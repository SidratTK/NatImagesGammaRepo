% generates Fig 5. 
% typical electrode for one monkey
% displays image patches at Receptive field and approximated patches
% also displays actual LFP power vs frequency, and 
% actual vs estimated gamma band response

function getTypicalElecFig

%1. fix some variables to be used later
folderSourceString= pwd; 
posList       = 1;    % which protocol? 1 for aH Flr
fValsToUse    = 17:32;% image numbers 
channelNumber = 34;   % elec to plot 
powerOption   = 2;

patchSizeDeg    = 3;            % images cropped to square patches, halflength
radiusMatrixDeg = 0.3:0.3:2.1;  % these sizes will be checked

selectOptions.meanThr =  [0.05 0.1 0.1]; 
selectOptions.stdThr  = 2*selectOptions.meanThr;
selectOptions.measure = 'diff';
selectOptions.method  = 'vector';

% 2. load details
[experimentalDetails,~] = getExperimentalDetails;
tmp = experimentalDetails{posList(1)};
subjectName     = tmp{1};
imageFolderName = tmp{4};
expDate         = tmp{5};
protocolName    = tmp{6};
rawImageFolder  = fullfile(folderSourceString,'data','images',imageFolderName);
numStimuli = length(fValsToUse);

% 3. make figure panels
G = figure; colormap gray
G.Units = 'centimeters'; G.Name = [subjectName,' elec',num2str(channelNumber),' poweroption ', num2str(powerOption)];
G.PaperType = 'a4';    G.PaperUnits = 'centimeters'; 
G.PaperSize = [17 9];  G.PaperOrientation = 'Portrait'; 
G.PaperPosition = [0 0 G.PaperSize]; 
G.Color = [1 1 1];     G.Position = [0 0 G.PaperSize]; 

hImagesPlot       = getPlotHandles(1,numStimuli,[0.04 0.83 0.75 0.1],0.002);
hImagePatchesPlot = getPlotHandles(1,numStimuli,[0.04 0.67 0.75 0.1],0.002);
hImagePatchPredictionPlot = getPlotHandles(1,numStimuli,[0.04 0.51 0.75 0.1],0.002);
hDataPlot         = getPlotHandles(1,numStimuli,[0.04 0.12 0.75 0.2],0.002);
hPowerPredictionPlot = subplot('Position',[0.85 0.58 0.14 0.30]);
hCorrelationPlot  = subplot('Position',[0.85 0.12 0.14 0.25]);
plotColor = 'g'; colorNamesPower = jet(numStimuli);
bColor = 'br'; bWidth = [0.9 0.6]; tickLength=[0.05 0.05];
fontSizeTicks = 9; fontSizeSmall = 10; fontSizeLarge = 11;

% 4. Get actual power
[powerST,powerBL,electrodeListPower,~,psdST,psdBL,freqVals] = getMeanEnergy(subjectName,expDate,protocolName,folderSourceString);
if powerOption==1
    % Do nothing. PowerST will be used.
elseif powerOption==2
    powerST = powerST ./ powerBL; % Ratio between ST and BL
elseif powerOption==3
    [powerST2,powerBL2] = getMeanEnergy(subjectName,expDate,protocolName,folderSourceString,{[80 150]}); % Take power between 80 to 150 Hz
    powerST = powerST ./ powerBL; % Ratio between ST and BL
    powerST2 = powerST2 ./ powerBL2; % Ratio between ST and BL in high gamma
    powerST = powerST - powerST2; % Difference in the ratios
end
% 4b. show LFP psd plots
psdST = squeeze(psdST(electrodeListPower==channelNumber,fValsToUse,:));
psdBL = squeeze(psdBL(electrodeListPower==channelNumber,fValsToUse,:));
for im=1:numStimuli
    hold(hDataPlot(im),'on');
    xline(hDataPlot(im),[30 80 150],'color',[0.7 0.7 0.7]); % to get the ranges of gamma & hi gamma 
    plot(hDataPlot(im),freqVals,zeros(1,length(freqVals)),'color',[0.5 0.5 0.5]);
    plot(hDataPlot(im),freqVals,log10(psdST(im,:)./psdBL(im,:)),'color',colorNamesPower(im,:));
    xlim(hDataPlot(im),[0 155]);
    ylim(hDataPlot(im),[-0.5 1.1]);
    set(hDataPlot(im),'XTick',[0 150],'XTickLabel',[],'YTick',[0 1],'YTickLabel','','tickdir','in','tickLength',tickLength,'FontSize',fontSizeTicks,'FontWeight','bold');
end 
set(hDataPlot(1),'XTick',[0 150],'XTickLabel',[0 150],'YTick',[0 1],'YTickLabel',[0 1],'XTickLabelRotation',0);
text(-0.3,-0.4,'Frequency','FontSize',fontSizeTicks,'FontWeight','bold','Parent',hDataPlot(1),'units','normalized');
text(-0.65,0,'Log Power','Rotation',90,'FontSize',fontSizeTicks,'FontWeight','bold','Parent',hDataPlot(1),'units','normalized');

% 5. get image patches and parameters
if isempty(radiusMatrixDeg)
    radiusMatrixDeg = 0.3:0.3:patchSizeDeg;
end
allStimParams = plotImageData2(hImagesPlot,hImagePatchesPlot,hImagePatchPredictionPlot,rawImageFolder,fValsToUse,channelNumber,subjectName,plotColor,selectOptions,patchSizeDeg,radiusMatrixDeg);
text(-0.3,-0.6,'Degrees','FontSize',fontSizeTicks,'FontWeight','bold','Parent',hImagePatchPredictionPlot(1),'units','normalized');
set(hImagePatchesPlot(1),'XTick',[-patchSizeDeg 0 patchSizeDeg],'XTickLabel',[-patchSizeDeg 0 patchSizeDeg],'YTick',[-patchSizeDeg 0 patchSizeDeg],'YTickLabel',[-patchSizeDeg 0 patchSizeDeg],'fontSize',fontSizeTicks,'FontWeight','bold')
set(hImagePatchPredictionPlot(1),'XTick',[-patchSizeDeg 0 patchSizeDeg],'XTickLabel',[-patchSizeDeg 0 patchSizeDeg],'YTick',[-patchSizeDeg 0 patchSizeDeg],'YTickLabel',[-patchSizeDeg 0 patchSizeDeg],'fontSize',fontSizeTicks,'FontWeight','bold')

% 6. get correlations 
allPower = squeeze(powerST(:,electrodeListPower==channelNumber,fValsToUse)); % power
[correlationsFull,correlationsSelected,predictionString,predictedPower,selectedImageIndices] = getAllCorrelations(subjectName,allStimParams,allPower);
hold(hPowerPredictionPlot,'on');
for i=1:numStimuli
    title(hImagesPlot(i),num2str(i),'color',colorNamesPower(i,:));
    if isempty(intersect(i,selectedImageIndices)) % Not a selected image
        plot(hPowerPredictionPlot,allPower(i),predictedPower(i),'marker','o','color',colorNamesPower(i,:));
    else 
        plot(hPowerPredictionPlot,allPower(i),predictedPower(i),'marker','o','color',colorNamesPower(i,:),'markerfacecolor',colorNamesPower(i,:));
    end 
end 
set(hPowerPredictionPlot,'YTickLabelRotation',90,'FontSize',fontSizeTicks,'FontWeight','bold','tickdir','in','tickLength',tickLength);
text(0.0,1.3,['el ' num2str(channelNumber)],'FontSize',fontSizeSmall,'FontWeight','bold','Parent',hPowerPredictionPlot,'units','normalized');
text(0.0,1.2,['rFull: ' num2str(round(correlationsFull(end),2))],'Color',bColor(1),'FontSize',fontSizeSmall,'FontWeight','bold','Parent',hPowerPredictionPlot,'units','normalized');
text(0.0,1.1,['rSel(N=' num2str(length(selectedImageIndices)) '):' num2str(round(correlationsSelected(end),2))],'Color',bColor(2),'FontSize',fontSizeSmall,'FontWeight','bold','Parent',hPowerPredictionPlot,'units','normalized');
text(-0.25,-0.2,'Actual Gamma Ratio','FontSize',fontSizeTicks,'FontWeight','bold','Parent',hPowerPredictionPlot,'units','normalized'); 
text(-0.25,-0.1,'Scaled Prediction','Rotation',90,'FontSize',fontSizeTicks,'FontWeight','bold','Parent',hPowerPredictionPlot,'units','normalized'); 
hold(hCorrelationPlot,'on');
bar(correlationsFull,'FaceColor',bColor(1),'BarWidth',bWidth(1),'FaceAlpha',0.7,'EdgeColor','none','Parent',hCorrelationPlot);
bar(correlationsSelected,'FaceColor',bColor(2),'BarWidth',bWidth(2),'FaceAlpha',0.6,'EdgeColor','none','Parent',hCorrelationPlot);
set(hCorrelationPlot,'XTick',1:length(predictionString),'XTickLabel',predictionString,'XTickLabelRotation',90,'YTickLabelRotation',90,'FontSize',fontSizeTicks,'FontWeight','bold','tickLength',tickLength,'tickdir','in');
text(0.1,1.2,'Full Set','Color',bColor(1),'FontSize',fontSizeTicks,'FontWeight','bold','Parent',hCorrelationPlot,'units','normalized');
text(0.1,1.1,'Selected Set','Color',bColor(2),'FontSize',fontSizeTicks,'FontWeight','bold','Parent',hCorrelationPlot,'units','normalized');
text(-0.25,0.1,'Correlation','Rotation',90,'FontSize',fontSizeTicks,'FontWeight','bold','Parent',hCorrelationPlot,'units','normalized');
ylim(hCorrelationPlot,[-0.5 1]);
        
text(-0.4,1.4,'A','FontSize',fontSizeLarge,'FontWeight','bold','Parent',hImagesPlot(1),'units','normalized');
text(-0.4,1.2,'B','FontSize',fontSizeLarge,'FontWeight','bold','Parent',hDataPlot(1),'units','normalized');
text(-0.3,1.3,'C','FontSize',fontSizeLarge,'FontWeight','bold','Parent',hPowerPredictionPlot,'units','normalized');
text(-0.3,1.3,'D','FontSize',fontSizeLarge,'FontWeight','bold','Parent',hCorrelationPlot,'units','normalized');

end

function allStimParams = plotImageData2(hImagesPlot,hImagePatches,hImagePatchPredictionPlot,rawImageFolder,fValsToUse,channelNumber,subjectName,colorName,selectOptions,patchSizeDeg,radiusMatrixDeg)

plottingDetails.displayPlotsFlag=1;
plottingDetailsG.displayPlotsFlag=0;
patchSizeDegG = 7; % cropped size for grating pipeline. bigger to account for sigma drop off

% Setting up standard gaborStimulus which is a patch (spatial frequency of 0 and phase of 90 degrees)
gaborStim.orientationDeg=0; % Does not matter since SF is zero
gaborStim.spatialFreqCPD=0; % For color patch
gaborStim.azimuthDeg=0;
gaborStim.elevationDeg=0;
gaborStim.sigmaDeg=100000; % The program makeGaborStimulus actually produces Gabors. However, when sigma is extremely large, it is essentially a grating whose radius is defined by another parameter

numImages = length(fValsToUse);
allStimParams = cell(1,numImages);
foldersSourceString=pwd;
for i=1:numImages
    imageFileName = fullfile(rawImageFolder,['Image' num2str(fValsToUse(i)) '.png']);
    plottingDetails.hImagePlot=hImagesPlot(i);
    plottingDetails.hImagePatches=hImagePatches(i);
    plottingDetails.colorNames=colorName;
    [patchData,imageAxesDeg] = getImagePatches(imageFileName,channelNumber,subjectName,foldersSourceString,patchSizeDeg,plottingDetails);
    stimParams = getSingleImageParameters(rgb2hsv(patchData{1}),imageAxesDeg,[0 0],radiusMatrixDeg,selectOptions,0);

    [patchDataG,imageAxesDegG] = getImagePatches(imageFileName,channelNumber,subjectName,foldersSourceString,patchSizeDegG,plottingDetailsG);
    tmpGParams = getSingleImageParametersGrating(rgb2hsv(patchDataG{1}),imageAxesDegG,[0 0],radiusMatrixDeg);
    if tmpGParams.categoryGabor 
        disp(['Img',num2str(i),' SF:',num2str(tmpGParams.spatialFreqCPD),' Or:',num2str(tmpGParams.orientationDeg),' Sz:',num2str(tmpGParams.sigmaDeg),...
            ' C:',num2str(tmpGParams.contrastPC),' Ov:',num2str(tmpGParams.oriVar)]);
    end
    stimParams.gaborParams = tmpGParams;
    allStimParams{i} = stimParams;
    
    tmpGaborStim = gaborStim;
    if tmpGParams.categoryGabor
        tmpGaborStim.sigmaDeg      = stimParams.gaborParams.sigmaDeg;
        tmpGaborStim.spatialFreqCPD= stimParams.gaborParams.spatialFreqCPD;
        tmpGaborStim.orientationDeg= stimParams.gaborParams.orientationDeg;
        tmpGaborStim.contrastPC    = stimParams.gaborParams.contrastPC;
        tmpGaborStim.radiusDeg     = 3*stimParams.gaborParams.sigmaDeg;
        tmpGaborStim.spatialFreqPhaseDeg = stimParams.gaborParams.spatialFreqPhaseDeg;
        tmpGaborPatch = makeGaborStimulus(tmpGaborStim,imageAxesDeg.xAxisDeg,imageAxesDeg.yAxisDeg,0);
        imagesc([-patchSizeDeg patchSizeDeg],[-patchSizeDeg patchSizeDeg],tmpGaborPatch,'Parent',hImagePatchPredictionPlot(i));
        clim(hImagePatchPredictionPlot(i),[0 1]);
    else 
        tmpGaborStim.hueDeg = stimParams.hueDeg;
        tmpGaborStim.sat = stimParams.sat;
        tmpGaborStim.contrastPC = stimParams.contrastPC;
        tmpGaborStim.spatialFreqPhaseDeg = stimParams.spatialFreqPhaseDeg;
        tmpGaborStim.radiusDeg = stimParams.radiusDeg;
        tmpGaborPatch = makeGaborStimulus(tmpGaborStim,imageAxesDeg.xAxisDeg,imageAxesDeg.yAxisDeg,0);
        imagesc([-patchSizeDeg patchSizeDeg],[-patchSizeDeg patchSizeDeg],tmpGaborPatch,'Parent',hImagePatchPredictionPlot(i));
    end  
    set(hImagesPlot(i),'XTick',[],'YTick',[],'XTicklabel',[],'YTicklabel',[],'box','off');
    set(hImagePatches(i),'XTick',[],'YTick',[],'XTicklabel',[],'YTicklabel',[],'box','off');
    set(hImagePatchPredictionPlot(i),'XTick',[],'YTick',[],'XTicklabel',[],'YTicklabel',[],'box','off');
    
end
end

