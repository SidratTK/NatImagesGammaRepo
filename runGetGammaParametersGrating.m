% run getGammaParameters.m
% runs code to get gamma parameters for grating stimuli.
% uses those parameters to predict gamma to defined stimuli
% also makes figure that shows model fits between 2 folds, and correlation.

function [elCorrM,modelParamsSubjects,modelParamsFold,G] = runGetGammaParametersGrating(showFig,useMediansFlag)

if ~exist('showFig','var'), showFig = 1; end
if ~exist('useMediansFlag','var'), useMediansFlag=1; end% use median params for prediction

tBL=[-0.25 0];tST=[0.25 0.50];     
fBandGamma=[30 80];
fBad = [];  

subjectNames= {'alpaH','kesariH'};
HueFlag=0;

% make figure panels:
G = figure; 
G.Units = 'centimeters'; G.PaperType = 'a4'; G.PaperUnits = 'centimeters'; G.PaperSize = [17.5 10];
G.PaperOrientation = 'Portrait'; G.PaperPosition = [0 0 G.PaperSize]; G.Color = [1 1 1]; G.Position = [0 0 G.PaperSize]; 
ePlot(1:2,:,1) = getPlotHandles(2,3,[0.10 0.12 0.60 0.82], 0.04, 0.12);
rPlot = getPlotHandles(2,1,[0.83 0.12 0.15 0.84], 0.05, 0.10);
H = figure;
H.Units = 'centimeters'; H.PaperType = 'a4'; H.PaperUnits = 'centimeters'; H.PaperSize = [11.5 16];
H.PaperOrientation = 'Portrait'; H.PaperPosition = [0 0 H.PaperSize]; H.Color = [1 1 1]; H.Position = [0 0 H.PaperSize];
histPlot=cell(1,2);
for m=1:2  % for 2 monkeys
    histPlot{m} = getPlotHandles(4,2,[0.10+0.5*(m-1) 0.09 0.38 0.86],0.05,0.11);
    histPlot{m} = histPlot{m}';histPlot{m} = histPlot{m}(:);
end 
F = figure;
F.Units = 'centimeters'; F.PaperType = 'a4'; F.PaperUnits = 'centimeters'; F.PaperSize = [17 11];
F.PaperOrientation = 'Portrait'; F.PaperPosition = [0 0 F.PaperSize]; F.Color = [1 1 1]; F.Position = [0 0 F.PaperSize]; 
ePlot(1:2,:,2) = getPlotHandles(2,3,[0.07 0.12 0.50 0.82], 0.04, 0.12);

fontSizeLabel = 12; fontSizeTicks = 9; fontSizeSmall = 10; fontSizeLarge = 11;
ticklength = [0.03,0.1];
dd = 0.09; 
col = parula(4); 
elecCol   = [0.4 0.4 0.4; 1 0 0]; % 'kr' 
xlabs  ={'SF','Size','Contrast'};     
labs   = {'SF-Ori ','Size-Ori ','Con-Ori '};
mAlpha = [0.2,0.5];

% for 2 subjects
modelParamsSubjects =cell(1,2);  % has median parameters
modelParamsAll      =cell(1,2);  % has all parameters. for use in histogram
modelParamsFold     =cell(2,3);  % fold wise parameters. for use in prediction. 3 iterations
gammaVals           =cell(2,3);  % gamma power responses
gammaEst            =cell(2,3);  % predict gamma responses
elCorr              =cell(2,3,3);
elCorrSig           =cell(2,3,3);
elCorrM             =cell(2,3);

for m=1:2  
    subjectName = subjectNames{m};
    % make model:
    for it=1:3
        [modelParamsSubjects{m},modelParamsAll{m},modelParamsFold{m,it},gammaVals{m,it}]= getGammaParameters(subjectName,tBL,tST,fBandGamma,fBad,HueFlag);
    end
end
for m=1:2  
    subjectName = subjectNames{m};
     
    % show parameters in a histogram:
    for pp=1:8
        if     pp<=4,         elecs = gammaVals{m,1}.sfOri.electrodes; 
        elseif pp>4 && pp<=6, elecs = gammaVals{m,1}.sizeOri.electrodes;
        elseif pp>6 && pp<=8, elecs = gammaVals{m,1}.conOri.electrodes;
        end
        parm = modelParamsAll{m}{pp}; parmOther = modelParamsAll{mod(m,2)+1}{pp};
        n = linspace(min([parm(:);parmOther(:)])-0.1,max([parm(:);parmOther(:)])+0.1,25);
        parmL= parm(elecs<=81 ); 
        parmE= parm(elecs>81  ); 
        mednP= median([parmL(:);parmE(:)]);
        if pp==3  % if Ori center angle
            parmL(parmL<0)= 360+parmL(parmL<0); parmE(parmE<0)= 360+parmE(parmE<0);
            mednP = rad2deg(circ_median(deg2rad([parmL(:);parmE(:)])));
            mednP(mednP<0)= 360+mednP(mednP<0);
        end
        [hl,n]= histcounts(parmL,n); hold(histPlot{m}(pp),'on');
        bar(histPlot{m}(pp),(n(1:end-1)+diff(n)/2),hl,'FaceColor',elecCol(1,:),'FaceAlpha',0.7,'EdgeColor','none','BarWidth',0.9);
        [he,~]= histcounts(parmE,n);
        bar(histPlot{m}(pp),(n(1:end-1)+diff(n)/2),he,'FaceColor',elecCol(2,:),'FaceAlpha',0.7,'EdgeColor','none','BarWidth',0.9);
        scatter(histPlot{m}(pp),mednP, max(hl)+3,22,'v','filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','k');
        xt1 = round(linspace(n(1),n(end),3),2); xtl = xt1;
        if ~isempty(find(diff(xt1)==0, 1)),  xt1 = n(1); xtl=xt1;  end % if not monotonic
        set(histPlot{m}(pp),'XTick',xt1,'XTickLabels',xtl,'YTick',0:15:60,'YTickLabels',0:15:60,'tickdir','in','ticklength',ticklength,'FontSize',10,'FontWeight','bold');
        ylim(histPlot{m}(pp),[0 max(hl)+10]); xlim(histPlot{m}(pp),[n(1) n(end)]+[n(1)-n(2) n(2)-n(1)]);
        xlabel(histPlot{m}(pp),modelParamsAll{m}{9}{pp});
    end
    ylabel(histPlot{m}(5),'No. of electrodes');
    text(1,1.2,['M',num2str(m)],'FontSize',fontSizeLarge,'FontWeight','bold','Parent',histPlot{m}(1),'units','normalized');
    
    % which are the LFPs & Ecogs for FS & size protocols?
    lfpInds{1} = gammaVals{m,1}.sfOri.electrodes<=81   ;   
    lfpInds{2} = gammaVals{m,1}.sizeOri.electrodes<=81 ; 
    lfpInds{3} = gammaVals{m,1}.conOri.electrodes<=81  ;  
    ecogInds{1}= gammaVals{m,1}.sfOri.electrodes>81    ;    
    ecogInds{2}= gammaVals{m,1}.sizeOri.electrodes>81  ;  
    ecogInds{3}= gammaVals{m,1}.conOri.electrodes>81   ;   
                
    % get estimated model gamma:
    for f= 1:2  % for 2 folds
        for it=1:3
            tempP = modelParamsFold{m,it}.median(f,:);  % this is median across all elecs- lfp+ecog
            [modelParamsMedian.sfCenter,modelParamsMedian.sfSigma,modelParamsMedian.oriCenter,modelParamsMedian.oriSpreadK,...
                modelParamsMedian.sizeSlope,modelParamsMedian.sizeAtHalfMax,modelParamsMedian.conSlope,modelParamsMedian.conAtHalfMax] = deal(tempP{:});
    
            % predict sfori
            stimParams.spatialFreqCPD = repmat(gammaVals{m,it}.sfOri.spatialFreqCPD',[1 length(gammaVals{m,it}.sfOri.orientationDeg)]);stimParams.spatialFreqCPD =stimParams.spatialFreqCPD(:);
            stimParams.orientationDeg = repmat(gammaVals{m,it}.sfOri.orientationDeg,[length(gammaVals{m,it}.sfOri.spatialFreqCPD) 1]); stimParams.orientationDeg=stimParams.orientationDeg(:);
            stimParams.radiusDeg      = gammaVals{m,it}.sfOri.radiusDeg;  % 
            stimParams.contrastPC     = gammaVals{m,it}.sfOri.contrastPC;
            elecs = gammaVals{m,it}.sfOri.electrodes;
            for el = 1:length(elecs)
                dataGamma = squeeze(gammaVals{m,it}.sfOri.delGammaFold{mod(f,2)+1}(el,:,:)); 
                if useMediansFlag
                    temp = getPredictedGamma(subjectName,stimParams,modelParamsMedian,dataGamma(:));
                else % if individual params
                    tempP = cat(2,cellfun(@(x) x(el), modelParamsFold{m,it}.all(f,1:4),'un',0),cellfun(@(x) x(gammaVals{m,it}.sizeOri.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,5:6),'un',0),cellfun(@(x) x(gammaVals{m,it}.conOri.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,7:8),'un',0));
                    [modelParamsInd.sfCenter,modelParamsInd.sfSigma,modelParamsInd.oriCenter,modelParamsInd.oriSpreadK,modelParamsInd.sizeSlope,modelParamsInd.sizeAtHalfMax,modelParamsInd.conSlope,modelParamsInd.conAtHalfMax] = deal(tempP{:});
                    if isempty(modelParamsInd.sizeSlope), modelParamsInd.sizeSlope = modelParamsMedian.sizeSlope; modelParamsInd.sizeAtHalfMax = modelParamsMedian.sizeAtHalfMax; end
                    temp = getPredictedGamma(subjectName,stimParams,modelParamsInd,dataGamma(:));
                end
                gammaEst{m,it}.sfOri.delGammaFold{f}(el,:,:)= reshape(temp,[length(gammaVals{m,it}.sfOri.spatialFreqCPD) length(gammaVals{m,it}.sfOri.orientationDeg)]);
            end  
    
            % predict sizeori
            stimParams.orientationDeg= repmat(gammaVals{m,it}.sizeOri.orientationDeg,[length(gammaVals{m,it}.sizeOri.radiusDeg) 1]); stimParams.orientationDeg=stimParams.orientationDeg(:);
            stimParams.spatialFreqCPD= gammaVals{m,it}.sizeOri.spatialFreqCPD;
            stimParams.radiusDeg     = repmat(gammaVals{m,it}.sizeOri.radiusDeg',[1 length(gammaVals{m,it}.sizeOri.orientationDeg)]); stimParams.radiusDeg=stimParams.radiusDeg(:);
            stimParams.contrastPC    = gammaVals{m,it}.sizeOri.contrastPC;
            elecs = gammaVals{m,it}.sizeOri.electrodes;
            for el = 1:length(elecs)
                dataGamma = squeeze(gammaVals{m}.sizeOri.delGammaFold{mod(f,2)+1}(el,:,:)); 
                if useMediansFlag
                    temp = getPredictedGamma(subjectName,stimParams,modelParamsMedian,dataGamma(:));
                else  % if individual params
                    tempP = cat(2,cellfun(@(x) x(gammaVals{m,it}.sfOri.electrodes==elecs(el)), modelParamsFold{m,it}.all(f,1:4),'un',0),cellfun(@(x) x(el),modelParamsFold{m,it}.all(f,5:6),'un',0),cellfun(@(x) x(gammaVals{m,it}.conOri.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,7:8),'un',0));
                    [modelParamsInd.sfCenter,modelParamsInd.sfSigma,modelParamsInd.oriCenter,modelParamsInd.oriSpreadK,modelParamsInd.sizeSlope,modelParamsInd.sizeAtHalfMax,modelParamsInd.conSlope,modelParamsInd.conAtHalfMax] = deal(tempP{:});
                    temp = getPredictedGamma(subjectName,stimParams,modelParamsInd,dataGamma(:));
                end
                gammaEst{m,it}.sizeOri.delGammaFold{f}(el,:,:) = reshape(temp,[length(gammaVals{m,it}.sizeOri.radiusDeg) length(gammaVals{m,it}.sizeOri.orientationDeg)]);
            end
    
            % predict conori
            stimParams.orientationDeg = repmat(gammaVals{m,it}.conOri.orientationDeg,[length(gammaVals{m,it}.conOri.contrastPC) 1]);stimParams.orientationDeg=stimParams.orientationDeg(:);
            stimParams.radiusDeg      = gammaVals{m,it}.conOri.radiusDeg;  % 
            stimParams.spatialFreqCPD = gammaVals{m,it}.conOri.spatialFreqCPD;
            stimParams.contrastPC     = repmat(gammaVals{m,it}.conOri.contrastPC',[1 length(gammaVals{m,it}.conOri.orientationDeg)]); stimParams.contrastPC=stimParams.contrastPC(:);
            elecs = gammaVals{m,it}.conOri.electrodes;
            for el = 1:length(elecs)
                dataGamma = squeeze(gammaVals{m,it}.conOri.delGammaFold{mod(f,2)+1}(el,:,:)); 
                if useMediansFlag
                    temp      = getPredictedGamma(subjectName,stimParams,modelParamsMedian,dataGamma(:));
                else % if individual params
                    tempP = cat(2,cellfun(@(x) x(gammaVals{m,it}.sfOri.electrodes==elecs(el)), modelParamsFold{m,it}.all(f,1:4),'un',0),cellfun(@(x) x(gammaVals{m,it}.sizeOri.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,5:6),'un',0),cellfun(@(x) x(el),modelParamsFold{m,it}.all(f,7:8),'un',0));
                    [modelParamsInd.sfCenter,modelParamsInd.sfSigma,modelParamsInd.oriCenter,modelParamsInd.oriSpreadK,modelParamsInd.sizeSlope,modelParamsInd.sizeAtHalfMax,modelParamsInd.conSlope,modelParamsInd.conAtHalfMax] = deal(tempP{:});
                    if isempty(modelParamsInd.sizeSlope), modelParamsInd.sizeSlope = modelParamsMedian.sizeSlope; modelParamsInd.sizeAtHalfMax = modelParamsMedian.sizeAtHalfMax; end
                    temp = getPredictedGamma(subjectName,stimParams,modelParamsInd,dataGamma(:));
                end
                gammaEst{m,it}.conOri.delGammaFold{f}(el,:,:) = reshape(temp,[length(gammaVals{m,it}.conOri.contrastPC) length(gammaVals{m,it}.conOri.orientationDeg)]);
            end
        end 
    end

    % show data & estimate in figure
    for pp=1:3
        for it=1:3
            if pp==1
                useData = gammaVals{m,it}.sfOri;
                useEst  = mean(cat(4,gammaEst{m,it}.sfOri.delGammaFold{1},gammaEst{m,it}.sfOri.delGammaFold{2}),4); % mean over folds
                useXs   = useData.spatialFreqCPD; useXinds=1:length(useData.spatialFreqCPD);
                useXlabs= useData.spatialFreqCPD;
            elseif pp==2
                useData = gammaVals{m,it}.sizeOri;
                useEst  = mean(cat(4,gammaEst{m,it}.sizeOri.delGammaFold{1},gammaEst{m,it}.sizeOri.delGammaFold{2}),4);
                useXs   = 1:length(useData.radiusDeg); useXinds= 1:length(useData.radiusDeg);%unique([1:2:length(useData.radiusDeg) length(useData.radiusDeg)]); 
                useXlabs= useData.radiusDeg;
            elseif pp==3
                useData = gammaVals{m,it}.conOri;
                useEst  = mean(cat(4,gammaEst{m,it}.conOri.delGammaFold{1},gammaEst{m,it}.conOri.delGammaFold{2}),4);
                useXs   = useData.contrastPC; useXinds= 1:length(useData.contrastPC);%unique([1:2:length(useData.contrastPC) length(useData.contrastPC)]);
                useXlabs= useData.contrastPC;
            end
            % electrode wise correlation
            for el = 1:size(useEst,1)
                [elCorr{m,pp,it}(el),elCorrSig{m,pp,it}(el)] = corr(squeeze(useData.delGamma(el,:))',squeeze(useEst(el,:))');
            end
        end 
        elCorrM{m,pp} = mean(cell2mat(elCorr(m,pp,:)),3);
        
        useOrs  = useData.orientationDeg;
        useCols = hsv(length(useData.orientationDeg));
        useElecInds = {lfpInds{pp} , ecogInds{pp}};  % Lfp and ECoGs
        hold(rPlot(m),'on'); ym=0;
        for ee=1:2 
            hold(ePlot(m,pp,ee),'on');
            % show data and fits
            pltData = squeeze(mean(useData.delGamma(useElecInds{ee},:,:),1));  % mean across elecs
            pltEst  = squeeze(mean(useEst(useElecInds{ee},:,:),1));         
            hD = plot(ePlot(m,pp,ee),1:length(useXs),squeeze(pltData),'marker','o','linestyle','none');
            hE = plot(ePlot(m,pp,ee),1:length(useXs),pltEst,'LineWidth',2,'Color',[0.5 0.5 0.5],'marker','.');
            set(hD,{'color'},num2cell(useCols,2)); set(hE,{'color'},num2cell(useCols,2));
            set(ePlot(m,pp,ee),'XTick',(useXinds),'XTickLabels',round(useXlabs(useXinds),2),'XTickLabelRotation',45,'YTick',0:5:20,'YTickLabels',0:5:20,'tickdir','in','ticklength',ticklength,'FontWeight','bold','FontSize',fontSizeTicks); 
            ym= max([ym max([pltData(:);pltEst(:)])]);
            xlim(ePlot(m,pp,ee),[0.5 length(useXs)+1]);
            text(0.05,0.85,[' n=',num2str(sum(useElecInds{ee}))],'FontWeight','bold','FontSize',fontSizeTicks,'Color',elecCol(ee,:),'units','normalized','Parent',ePlot(m,pp,ee));
            if m==1
                title(ePlot(m,pp,ee),labs{pp},'FontWeight','bold','FontSize',fontSizeLarge);
            else
                xlabel(ePlot(m,pp,ee),xlabs{pp});
            end
            % show corr across elecs 
            useCorr = elCorrM{m,pp}(useElecInds{ee});
            scatter(rPlot(m),repmat(pp,size(useCorr)),useCorr,26,elecCol(ee,:),'filled','Marker','o','MarkerFaceAlpha',mAlpha(ee),'MarkerEdgeAlpha',1);
        end  
        useCorr = elCorrM{m,pp}; % all elecs
        mn = mean(useCorr); st = std(useCorr)./sqrt(length(useCorr));
        bar(rPlot(m),pp,mn,'EdgeColor',col(pp,:),'EdgeAlpha',1,'FaceAlpha',0,'LineWidth',2);
        errorbar(rPlot(m),pp,mn,st,'Color','k');
        text(0.8,1.1-pp*dd,['n=',num2str(length(useCorr))],'FontWeight','bold','FontSize',fontSizeTicks,'Color',col(pp,:),'units','normalized','Parent',rPlot(m));
    end
    ylim(rPlot(m),[0 1.3]); 
    ylabel(rPlot(m),'Correlation','FontSize',fontSizeSmall,'FontWeight','bold');
    text(-0.50,0.4, ['M ',num2str(m)],'Rotation',90,'FontWeight','bold','FontSize',fontSizeSmall,'units','normalized','Parent',rPlot(m));
    set(rPlot(m),'XTick',1:3,'XTickLabel',labs,'YTick',0:0.5:1,'YTickLabel',0:0.5:1,'tickdir','in','FontWeight','bold','FontSize',fontSizeTicks)
    xlim(rPlot(m),[0.5 3.5]);
    for oo=1:length(useOrs)
        text(0.92,1.0-dd*oo, num2str(useOrs(oo)),'FontWeight','bold','FontSize',fontSizeTicks,'units','normalized','Parent',ePlot(m,3,1),'Color',useCols(oo,:));
    end 
    text(-0.40,0.4, ['M ',num2str(m)],'Rotation',90,'FontWeight','bold','FontSize',fontSizeSmall,'units','normalized','Parent',ePlot(m,1,1));
    text(-0.25,0.1,'Gamma ratio','Rotation',90,'FontSize',fontSizeSmall,'FontWeight','bold','Units','normalized','Parent',ePlot(m,1,1));
    text(0.05, 0.95, 'LFP','FontWeight','bold','FontSize',fontSizeSmall,'units','normalized','Parent',ePlot(m,1,1));
    text(0.05, 0.95, 'ECoG','FontWeight','bold','FontSize',fontSizeSmall,'units','normalized','Parent',ePlot(m,1,2));
    for pp=1:3
        for ee=1:2
            ylim(ePlot(m,pp,ee),[0 ym]);
        end
    end
end 
text(-0.45, 1.1, 'A','FontWeight','bold','FontSize',fontSizeLabel,'units','normalized','Parent',ePlot(1,1,1));
text(-0.45, 1.05,'B','FontWeight','bold','FontSize',fontSizeLabel,'units','normalized','Parent',rPlot(1));
set(rPlot(1),'XTickLabel',[]);

if ~showFig
close(F); close(G);
close(H);
end
