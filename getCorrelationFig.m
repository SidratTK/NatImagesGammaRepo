% analyzes Data on all selected sessions.
% plots bar graph of averaged correlations between predicted and actual
% image gamma
% generates Figure 6.


folderSourceString= pwd;
patchSizeDeg = 3; % for patch
patchSizeDegG= 7; % for grating
radiusMatrixDeg = 0.3:0.3:2.1; % size steps for calculations

selectOptions.meanThr = [0.05 0.1 0.1];
selectOptions.stdThr = 2*selectOptions.meanThr;
selectOptions.measure = 'diff';
selectOptions.method = 'vector';

powerOptions = [2 2 3]; % 1 - ST power, 2 - ST/BL ratio, 3 - ST/BL ratio minus ST/BL ratio in HG
powerInRanges= {[30 80],[80 150],[30 80]};
runForInds   = [1 2 3];

cutoffSelectedImages = 8;

[experimentalDetails,matchIndex] = getExperimentalDetails;
posList = [5 2 1 10 13 14]; % Index for which data needs to be analyzed
setList = {[1 2] 1 [1 2] [1 2] 1 [1 2]}; % Two image sets per experimental session except "human" dataset

G = figure; 
G.Units = 'centimeters'; 
G.PaperType = 'a4';   G.PaperUnits = 'centimeters'; 
G.PaperSize = [17.5 16];  G.PaperOrientation = 'Portrait';  % G.PaperSize = [9 16];
G.PaperPosition = [0 0 G.PaperSize]; 
G.Color = [1 1 1];    G.Position = [0 0 G.PaperSize]; 
bPlots(1,:) = getPlotHandles(1,4,[0.1 0.64 0.85 0.24],0.06,0.05);
bPlots(2,:) = getPlotHandles(1,4,[0.1 0.37 0.85 0.24],0.06,0.05);
bPlots(3,:) = getPlotHandles(1,4,[0.1 0.10 0.85 0.24],0.06,0.05);

bCols = 'bmg'; bWidths=[0.8 0.8 0.8]; col2s = {[0.4 0.4 0.4],[1 0 0]}; 
mCols = {[0.7 0.7 0.7],[0.4 0.4 0.4]}; mMarks = 'os';

fontSizeLabel = 12; fontSizeTicks = 9; fontSizeSmall = 10; fontSizeLarge = 11;

powerST = cell(1,length(runForInds));
correlationsFullAll = cell(1,length(runForInds));
correlationsSelectedAll= cell(1,length(runForInds));
numSelectedPatches = []; 
electrodeListAll=[];
subInd = [];
selectedGAll=[]; selectedGratings=[];
allStimParamsMain= cell(length(posList),2);

for i=1:length(posList)
    tmp = experimentalDetails{posList(i)};
    
    for j=setList{i} % Two image sets per experimental session except "human" dataset
        subjectName = tmp{1};
        imageFolderName = tmp{4};
        expDate = tmp{5};
        protocolName = tmp{6};
        rawImageFolder = fullfile(folderSourceString,'data','images',imageFolderName);
        if j==1
            dataType = tmp{2};
            imageIndices = 1:16;
        else
            dataType = tmp{3};
            imageIndices = 17:32;
        end
        
        if ~isempty(dataType)
            disp(['Working on ' subjectName expDate protocolName ', set: ' dataType]);
            % get power for r=123
            for r = runForInds
                powerInRange= powerInRanges(r);
                powerOption = powerOptions(r);
                [powerST{r},powerBL,electrodeList] = getMeanEnergy(subjectName,expDate,protocolName,folderSourceString,powerInRange);
                if powerOption==1
                    powerST{r} = squeeze(powerST{r}(:,:,imageIndices)); % Only take stimulus power
                elseif powerOption==2
                    powerST{r} = squeeze(powerST{r}(:,:,imageIndices)) ./ squeeze(powerBL(:,:,imageIndices)); % Ratio between ST and BL
                elseif powerOption==3
                    [powerST2,powerBL2,electrodeList] = getMeanEnergy(subjectName,expDate,protocolName,folderSourceString,{[80 150]}); % Take power between 80 to 150 Hz
                    powerST{r} = squeeze(powerST{r}(:,:,imageIndices)) ./ squeeze(powerBL(:,:,imageIndices)); % Ratio between ST and BL
                    powerST2 = squeeze(powerST2(:,:,imageIndices)) ./ squeeze(powerBL2(:,:,imageIndices)); % Ratio between ST and BL in high gamma
                    powerST{r} = powerST{r} - powerST2;
                end  
            end 

            % Get stimulus parameters for image patches
            disp('Getting stim params...');
            numElectrodes = size(powerST{1},1);
            numImages = size(powerST{1},2);
            allStimParams = cell(numImages,numElectrodes);
            selectedG = false(numImages,numElectrodes);
            for im=1:numImages
                % Load image
                plottingDetails.displayPlotsFlag=0;
                imageFileName = fullfile(rawImageFolder,['Image' num2str(imageIndices(im)) '.png']);
                [patchData,imageAxesDeg] = getImagePatches(imageFileName,electrodeList,subjectName,folderSourceString,patchSizeDeg,plottingDetails);
                [patchDataG,imageAxesDegG]=getImagePatches(imageFileName,electrodeList,subjectName,folderSourceString,patchSizeDegG,plottingDetails);
                % Get Stim Parameters
                for el=1:numElectrodes
                    stimParams = getSingleImageParameters(rgb2hsv(patchData{el}),imageAxesDeg,[0 0],radiusMatrixDeg,selectOptions,0);
                    stimGabParams = getSingleImageParametersGrating(rgb2hsv(patchDataG{el}),imageAxesDegG,[0 0],radiusMatrixDeg);
                    stimParams.gaborParams = stimGabParams;
                    allStimParams{im,el} = stimParams;
                    selectedG(im,el)=stimGabParams.categoryGabor;
                end
            end
            selectedGAll = cat(2,selectedGAll,selectedG);
            selectedGratings   = cat(2,selectedGratings,sum(selectedG,1));
            
            % get correlations
            disp('Getting correlations');
            for r = runForInds
                cFull = zeros(7,numElectrodes);
                cSelected = zeros(7,numElectrodes);
                numSI = zeros(1,numElectrodes);
                for el=1:numElectrodes
                    actualPower = powerST{r}(el,:);
                    stimParams  = allStimParams(:,el);
                    [cFull(:,el),cSelected(:,el),predictionString,~,selectedImageIndices] = getAllCorrelations(subjectName,stimParams,actualPower);
                    numSI(el) = length(selectedImageIndices);
                end
                correlationsFullAll{r}= cat(2,correlationsFullAll{r},cFull);
                correlationsSelectedAll{r} = cat(2,correlationsSelectedAll{r},cSelected);
            end
            numSelectedPatches= cat(2,numSelectedPatches,numSI);
            electrodeListAll = cat(1,electrodeListAll,electrodeList(:));
            if     strcmp(subjectName,'alpaH')
                   subInd = cat(2,subInd,ones(1,numElectrodes));
            elseif strcmp(subjectName,'kesariH')
                   subInd = cat(2,subInd,2*ones(1,numElectrodes));
            end
        end 
    end
end

disp('Average number of selected grating patches are: ')
disp([num2str(mean(selectedGratings)),' +- sd ', num2str(std(selectedGratings))]);
disp('Average number of selected uniform patches are: ')
disp([num2str(mean(numSelectedPatches)),' +- sd', num2str(std(numSelectedPatches))]);

% plot bar graphs
for r = runForInds
    for m=1:2
        useInds= subInd==m;
        for i=1:2 % full and selected
            if i==1
                goodPos = useInds;
                X = correlationsFullAll{r}(:,goodPos);
                N = size(X,2);
                titleStr1 = 'Full set';
                titleStr2 = ['(N=' num2str(N) ')'];
            else 
                goodPos = useInds & (numSelectedPatches>=cutoffSelectedImages) ;
                X = correlationsSelectedAll{r}(:,goodPos);
                N = size(X,2);
                titleStr1 = 'Selected set';
                titleStr2 = ['(cutoff=' num2str(cutoffSelectedImages) ', N=' num2str(N) ')'];
            end 
            mData = mean(X,2);
            sData = std(X,[],2)/sqrt(N);
            ecogInd=electrodeListAll(goodPos)>81;
            hPlot = bPlots(r,i+m*(m-1)); hold(hPlot,'on');
            scatter(hPlot,1:length(mData),X(:,~ecogInd),14,'o','MarkerFaceColor',col2s{1},'MarkerEdgeColor',col2s{1},'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0);
            scatter(hPlot,1:length(mData),X(:,ecogInd),14,'o','MarkerFaceColor',col2s{2},'MarkerEdgeColor',col2s{2},'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0);
            bar(mData,'Parent',hPlot,'FaceColor',bCols(r),'FaceAlpha',0,'EdgeColor',bCols(r),'EdgeAlpha',1,'BarWidth',bWidths(r),'LineWidth',2);
            errorbar(hPlot,mData,sData,'k.');
            set(hPlot,'XTick',1:length(mData),'XTickLabel',[],'YTick',-0.5:0.5:1,'YTickLabel',-0.5:0.5:1,'fontSize',fontSizeTicks,'fontWeight','bold','box','off');
            ylim(hPlot,[-1 1]);
            if r==3
                set(hPlot,'XTick',1:length(mData),'XTickLabel',predictionString,'XTickLabelRotation',90);
            end
            if r==1
                text(0.1,1.20,titleStr1,'FontSize',fontSizeSmall,'FontWeight','bold','Parent',hPlot,'units','normalized');
                text(0.1,1.10,titleStr2,'FontSize',fontSizeSmall,'FontWeight','bold','Parent',hPlot,'units','normalized');
            end 
        end
    end
    text(-0.3,0.35,'Correlation','Rotation',90,'fontSize',fontSizeSmall,'fontWeight','bold','Parent',bPlots(r,1),'units','normalized');
end
text(-0.5,0.35,num2str(powerInRanges{1}),'Rotation',90,'color',bCols(1),'fontSize',fontSizeSmall,'fontWeight','bold','Parent',bPlots(1,1),'units','normalized')    
text(-0.5,0.35,num2str(powerInRanges{2}),'Rotation',90,'color',bCols(2),'fontSize',fontSizeSmall,'fontWeight','bold','Parent',bPlots(2,1),'units','normalized')    
text(-0.5,0.35,'Difference',             'Rotation',90,'color',bCols(3),'fontSize',fontSizeSmall,'fontWeight','bold','Parent',bPlots(3,1),'units','normalized')    
text(1.1,1.35,'M1','FontSize',fontSizeLarge,'FontWeight','bold','Parent',bPlots(1,1),'units','normalized');
text(1.1,1.35,'M2','FontSize',fontSizeLarge,'FontWeight','bold','Parent',bPlots(1,3),'units','normalized');
text(0.7,-0.3,'Parameters chosen','fontSize',fontSizeLarge,'fontWeight','bold','Parent',bPlots(3,1),'units','normalized');
text(0.7,-0.3,'Parameters chosen','fontSize',fontSizeLarge,'fontWeight','bold','Parent',bPlots(3,3),'units','normalized');


