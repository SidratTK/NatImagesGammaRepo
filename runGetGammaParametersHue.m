% run getGammaParameters.m
% gets gamma parameters for hue based model.
% uses these parameters to predict gamma response to defined stimuli
% also makes figure that shows model fits

function [elCorrM,modelParamsSubjects,modelParamsFold,figM] = runGetGammaParametersHue(showFig,useMediansFlag)

if ~exist('showFig','var'), showFig = 1; end
if ~exist('useMediansFlag','var'), useMediansFlag=1; end% use median params for prediction

tBL=[-0.25 0]; tST=[0.25 0.50];  
fBandGamma=[30 80];
fBad = [];  

HueFlag=1;
subjectNames= {'alpa','kesari'};

% make figure panels:
figM = figure;  % for data gamma
figM.Units      = 'centimeters'; figM.PaperType  = 'a4';
figM.PaperUnits = 'centimeters'; figM.PaperSize  = [17.5 11]; %
figM.PaperPosition = [0 0 figM.PaperSize]; figM.Color = [1 1 1];
figM.Position      = [0 0 figM.PaperSize]; figM.PaperOrientation = 'Portrait';
ePlot = getPlotHandles(2,4,[0.08 0.14 0.66 0.80],0.04,0.12);
rPlot = getPlotHandles(2,1,[0.83 0.14 0.15 0.80],0.04,0.10);
    
figH = figure;  % for histograms of parameters 
figH.Units      = 'centimeters'; figH.PaperType  = 'a4';
figH.PaperUnits = 'centimeters'; figH.PaperSize  = [17 17]; %
figH.PaperPosition = [0 0 figH.PaperSize]; figH.Color = [1 1 1];
figH.Position      = [0 0 figH.PaperSize]; figH.PaperOrientation = 'Portrait';
histPlot = cell(1,2);
for ii=1:2  % for 2 monkeys
        histPlot{ii} = getPlotHandles(5,3,[0.10+0.5*(ii-1) 0.07 0.38 0.89],0.06,0.09);
        histPlot{ii} = histPlot{ii}'; histPlot{ii} = histPlot{ii}(:); % rearrange
        set(histPlot{ii}([3,9,12,15]),'visible','off');
        histPlot{ii} = histPlot{ii}([1:2,4:8,10:11,13:14]);
end 
fontSizeLabel = 12; fontSizeTicks = 9; fontSizeSmall = 10; fontSizeLarge = 11;
col   = parula(5); dd=0.09; ticklength = [0.03,0.1];
hcols = hsv(360);
labs2 = {'Hue','Hue-Size','Hue-Sat','Hue-Val'};
xlabs = {'Hue','Size','Sat','Val'};

% for 2 subjects
modelParamsSubjects =cell(1,2);  % has median parameters
modelParamsAll      =cell(1,2);  % has all parameters. for use in histogram
modelParamsFold     =cell(2,3);  % fold wise parameters. for us ein prediction
gammaVals           =cell(2,3);  % gamma power responses
gammaEst            =cell(2,3);  % predict gamma responses
elCorr              =cell(2,4,3);
elCorrSig           =cell(2,4,3);
elCorrM             =cell(2,4);

for m=1:2
    subjectName = subjectNames{m};
    % make model:
    for it=1:3   % 3 iterations for folds
    [modelParamsSubjects{m},modelParamsAll{m},modelParamsFold{m,it},gammaVals{m,it}]= getGammaParameters(subjectName,tBL,tST,fBandGamma,fBad,HueFlag);
    end
end
for m=1:2
    subjectName = subjectNames{m};
        % show parameters in a histogram:
        for pp=1:11
            parm = modelParamsAll{m}{pp};
            mednP= median(parm);
            if pp==1||pp==3  % if Hue center angle
                parm(parm<0)= 2*pi+parm(parm<0);
                mednP = circ_median(parm(~isnan(parm)));
                mednP(mednP<0)= 2*pi+mednP(mednP<0);
            end 
            parmOther = modelParamsAll{mod(m,2)+1}{pp};
            n = linspace(0,max([parm(:);parmOther(:)])+0.1,25);
            [h,n] = histcounts(parm,n); hold(histPlot{m}(pp),'on');
            bar(histPlot{m}(pp),(n(1:end-1)+diff(n)/2),h,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.7,'EdgeColor','none','BarWidth',0.9);
            scatter(histPlot{m}(pp),mednP, max(h)+3,22,'v','filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','k');
            xt1=round(linspace(n(1),n(end),3),2); xtl = xt1;
            if ~isempty(find(diff(xt1)==0, 1)),  xt1 = n(1); xtl=xt1;  end % if not monotonic
            set(histPlot{m}(pp),'XTick',xt1,'XTickLabels',xtl,'YTick',0:15:60,'YTickLabels',0:15:60,'tickdir','in','ticklength',ticklength,'FontSize',fontSizeTicks,'FontWeight','bold');
            ylim(histPlot{m}(pp),[0 max(h)+10]); xlim(histPlot{m}(pp),[n(1) n(end)]+[n(1)-n(2) n(2)-n(1)]);
            xlabel(histPlot{m}(pp),modelParamsAll{m}{12}{pp});
        end
        ylabel(histPlot{m}(6),'No. of electrodes');
        text(1,1.2,['M',num2str(m),'R'],'FontSize',fontSizeLarge,'FontWeight','bold','Parent',histPlot{m}(1),'units','normalized');
    
    % get estimated model gamma:
    for f= 1:2  % for 2 folds
        for it=1:3
        tempP = modelParamsFold{m,it}.median(f,:);
        [modelParamsMedian.centerR,modelParamsMedian.spreadR,modelParamsMedian.centerC,modelParamsMedian.spreadC,modelParamsMedian.gainCbyR,modelParamsMedian.sizeSlope,...
            modelParamsMedian.sizeAtHalfMax,modelParamsMedian.satSlope,modelParamsMedian.satAtHalfMax,modelParamsMedian.valSlope,modelParamsMedian.valIntercept] = deal(tempP{:});

        stimParams.spatialFreqPhaseDeg = -90; % for use with hue patches
        % for hueScreen parameters
        stimParams.hueDeg    = gammaVals{m,it}.hueScreen.hueDeg;
        stimParams.radiusDeg = gammaVals{m,it}.hueScreen.radiusDeg;  % full screen
        stimParams.sat       = gammaVals{m,it}.hueScreen.sat;
        stimParams.contrastPC= gammaVals{m,it}.hueScreen.val * 100;
        elecs                = gammaVals{m,it}.hueScreen.electrodes;
        for el = 1:length(elecs)
            dataGamma = (gammaVals{m,it}.hueScreen.delGammaFold{mod(f,2)+1}(el,:)); 
            if useMediansFlag
                gammaEst{m,it}.hueScreen.delGammaFold{f}(el,:,:)= getPredictedGamma(subjectName,stimParams,modelParamsMedian,dataGamma);
            else
                tempP = cat(2,cellfun(@(x) x(el), modelParamsFold{m,it}.all(f,1:5),'un',0),cellfun(@(x) x(gammaVals{m,it}.hueSize.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,6:7),'un',0),...
                              cellfun(@(x) x(gammaVals{m,it}.hueSatFullscreen.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,8:9),'un',0),cellfun(@(x) x(gammaVals{m,it}.hueValFullscreen.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,10:11),'un',0));
                [modelParamsInd.centerR,modelParamsInd.spreadR,modelParamsInd.centerC,modelParamsInd.spreadC,modelParamsInd.gainCbyR,modelParamsInd.sizeSlope,...
                            modelParamsInd.sizeAtHalfMax,modelParamsInd.satSlope,modelParamsInd.satAtHalfMax,modelParamsInd.valSlope,modelParamsInd.valIntercept] = deal(tempP{:});
                if isempty(modelParamsInd.sizeSlope), modelParamsInd.sizeSlope = modelParamsMedian.sizeSlope; modelParamsInd.sizeAtHalfMax = modelParamsMedian.sizeAtHalfMax;  % if no size data, use median
                end 
                gammaEst{m,it}.hueScreen.delGammaFold{f}(el,:,:)= getPredictedGamma(subjectName,stimParams,modelParamsInd,dataGamma);
            end
        end 

        % for hueSize parameters
        stimParams.radiusDeg = gammaVals{m,it}.hueSize.radiusDeg'; stimParams.radiusDeg=stimParams.radiusDeg(:);
        stimParams.hueDeg    = zeros(length(gammaVals{m,it}.hueSize.radiusDeg),1); stimParams.hueDeg=stimParams.hueDeg(:); % 0deg = red
        stimParams.sat       = gammaVals{m,it}.hueSize.sat;
        stimParams.contrastPC= gammaVals{m,it}.hueSize.val * 100;
        elecs = gammaVals{m,it}.hueSize.electrodes;
        for el = 1:length(elecs) 
            dataGamma = squeeze(gammaVals{m,it}.hueSize.delGammaFold{mod(f,2)+1}(el,:,gammaVals{m,it}.hueSize.hueDeg==0)); 
            if useMediansFlag
                temp = getPredictedGamma(subjectName,stimParams,modelParamsMedian,dataGamma);
            else
                tempP = cat(2,cellfun(@(x) x(gammaVals{m,it}.hueScreen.electrodes==elecs(el)), modelParamsFold{m,it}.all(f,1:5),'un',0),cellfun(@(x) x(el),modelParamsFold{m,it}.all(f,6:7),'un',0),...
                              cellfun(@(x) x(gammaVals{m,it}.hueSatFullscreen.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,8:9),'un',0),cellfun(@(x) x(gammaVals{m,it}.hueValFullscreen.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,10:11),'un',0));
                [modelParamsInd.centerR,modelParamsInd.spreadR,modelParamsInd.centerC,modelParamsInd.spreadC,modelParamsInd.gainCbyR,modelParamsInd.sizeSlope,...
                            modelParamsInd.sizeAtHalfMax,modelParamsInd.satSlope,modelParamsInd.satAtHalfMax,modelParamsInd.valSlope,modelParamsInd.valIntercept] = deal(tempP{:});
                temp = getPredictedGamma(subjectName,stimParams,modelParamsInd,dataGamma);
            end
                gammaEst{m,it}.hueSize.delGammaFold{f}(el,:,:) = temp; % reshape(temp,[length(gammaVals.hueSize.radDeg) length(gammaVals.hueSize.hueDeg)]); 
        end

        % for Hue-Saturation
        stimParams.sat       = repmat(gammaVals{m,it}.hueSatFullscreen.sat,[1 length(gammaVals{m,it}.hueSatFullscreen.hueDeg)]); stimParams.sat=stimParams.sat(:);
        stimParams.hueDeg    = repmat(gammaVals{m,it}.hueSatFullscreen.hueDeg,[length(gammaVals{m,it}.hueSatFullscreen.sat) 1]); stimParams.hueDeg=stimParams.hueDeg(:);
        stimParams.radiusDeg = gammaVals{m,it}.hueSatFullscreen.radiusDeg;
        stimParams.contrastPC= gammaVals{m,it}.hueSatFullscreen.val * 100;
        elecs = gammaVals{m,it}.hueSatFullscreen.electrodes;
        for el = 1:length(elecs) 
            dataGamma = squeeze(gammaVals{m,it}.hueSatFullscreen.delGammaFold{mod(f,2)+1}(el,:,:)); 
            if useMediansFlag
                temp = getPredictedGamma(subjectName,stimParams,modelParamsMedian,dataGamma);
            else
                tempP = cat(2,cellfun(@(x) x(gammaVals{m,it}.hueScreen.electrodes==elecs(el)), modelParamsFold{m,it}.all(f,1:5),'un',0),cellfun(@(x) x(gammaVals{m,it}.hueSize.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,6:7),'un',0),...
                              cellfun(@(x) x(el),modelParamsFold{m,it}.all(f,8:9),'un',0),cellfun(@(x) x(gammaVals{m,it}.hueValFullscreen.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,10:11),'un',0));
                [modelParamsInd.centerR,modelParamsInd.spreadR,modelParamsInd.centerC,modelParamsInd.spreadC,modelParamsInd.gainCbyR,modelParamsInd.sizeSlope,...
                            modelParamsInd.sizeAtHalfMax,modelParamsInd.satSlope,modelParamsInd.satAtHalfMax,modelParamsInd.valSlope,modelParamsInd.valIntercept] = deal(tempP{:});
                if isempty(modelParamsInd.sizeSlope), modelParamsInd.sizeSlope = modelParamsMedian.sizeSlope; modelParamsInd.sizeAtHalfMax = modelParamsMedian.sizeAtHalfMax; 
                end
                temp = getPredictedGamma(subjectName,stimParams,modelParamsInd,dataGamma);
            end
            gammaEst{m,it}.hueSatFullscreen.delGammaFold{f}(el,:,:) = reshape(temp,[length(gammaVals{m,it}.hueSatFullscreen.sat) length(gammaVals{m,it}.hueSatFullscreen.hueDeg)]);
        end
        
        % for Hue Value stimuli
        stimParams.contrastPC = 100*repmat(gammaVals{m,it}.hueValFullscreen.val,[1 length(gammaVals{m,it}.hueValFullscreen.hueDeg)]);stimParams.contrastPC=stimParams.contrastPC(:);
        stimParams.hueDeg     = repmat(gammaVals{m,it}.hueValFullscreen.hueDeg,[length(gammaVals{m,it}.hueValFullscreen.val) 1]);    stimParams.hueDeg=stimParams.hueDeg(:);
        stimParams.sat        = gammaVals{m,it}.hueValFullscreen.sat; 
        stimParams.radiusDeg  = gammaVals{m,it}.hueValFullscreen.radiusDeg;
        elecs = gammaVals{m,it}.hueValFullscreen.electrodes;
        for el = 1:length(elecs)             
            dataGamma = squeeze(gammaVals{m,it}.hueValFullscreen.delGammaFold{mod(f,2)+1}(el,:,:)); 
            if useMediansFlag
                temp = getPredictedGamma(subjectName,stimParams,modelParamsMedian,dataGamma);
            else
                tempP = cat(2,cellfun(@(x) x(gammaVals{m,it}.hueScreen.electrodes==elecs(el)), modelParamsFold{m,it}.all(f,1:5),'un',0),cellfun(@(x) x(gammaVals{m,it}.hueSize.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,6:7),'un',0),...
                              cellfun(@(x) x(gammaVals{m,it}.hueSatFullscreen.electrodes==elecs(el)),modelParamsFold{m,it}.all(f,8:9),'un',0),cellfun(@(x) x(el),modelParamsFold{m,it}.all(f,10:11),'un',0));
                [modelParamsInd.centerR,modelParamsInd.spreadR,modelParamsInd.centerC,modelParamsInd.spreadC,modelParamsInd.gainCbyR,modelParamsInd.sizeSlope,...
                            modelParamsInd.sizeAtHalfMax,modelParamsInd.satSlope,modelParamsInd.satAtHalfMax,modelParamsInd.valSlope,modelParamsInd.valIntercept] = deal(tempP{:});
                if isempty(modelParamsInd.sizeSlope), modelParamsInd.sizeSlope = modelParamsMedian.sizeSlope; modelParamsInd.sizeAtHalfMax = modelParamsMedian.sizeAtHalfMax; end
                temp = getPredictedGamma(subjectName,stimParams,modelParamsInd,dataGamma);
            end
            gammaEst{m,it}.hueValFullscreen.delGammaFold{f}(el,:,:) = reshape(temp,[length(gammaVals{m,it}.hueValFullscreen.val) length(gammaVals{m,it}.hueValFullscreen.hueDeg)]);
        end
        end
    end

    % show data & estimate in a figure
    ym=0;
    for pp=[1 3 4 2]
        for it=1:3
            if pp==1
                useData = gammaVals{m,it}.hueScreen;
                useEst  = mean(cat(3,gammaEst{m,it}.hueScreen.delGammaFold{1},gammaEst{m,it}.hueScreen.delGammaFold{2}),3); % mean over folds
                useXs   = 1:length(useData.hueDeg); useXinds=1:6:36;
                useXtLabs= circshift(useData.hueDeg',size(useData.hueDeg,2)/2);
                hcols1  = hcols(intersect(0:length(hcols)-1,useData.hueDeg,'stable')+1,:);
                hcols1  = circshift(hcols1,size(hcols1,1)/2);
                pltData = squeeze(mean(useData.delGamma,1));  % mean across elecs
                pltData = circshift(pltData(:),size(pltData,2)/2);  % 
                pltEst  = squeeze(mean(useEst,1));   % mean across elecs
                pltEst  = circshift(pltEst(:),length(pltEst)/2);    %
            else  
                if pp==2
                    useData = gammaVals{m,it}.hueSize;
                    useEst  = mean(cat(4,gammaEst{m,it}.hueSize.delGammaFold{1},gammaEst{m,it}.hueSize.delGammaFold{2}),4);
                    useXs   = 1:length(useData.radiusDeg); useXinds= 1:length(useData.radiusDeg); % useXs = (useData.radiusDeg);
                    useXtLabs= useData.radiusDeg;
                elseif  pp==3
                    useData = gammaVals{m,it}.hueSatFullscreen;
                    useEst  = mean(cat(4,gammaEst{m,it}.hueSatFullscreen.delGammaFold{1},gammaEst{m,it}.hueSatFullscreen.delGammaFold{2}),4);
                    useXs   = 1:length(useData.sat); useXinds= 1:length(useData.sat);
                    useXtLabs= useData.sat;
                elseif pp==4
                    useData = gammaVals{m,it}.hueValFullscreen;
                    useEst  = mean(cat(4,gammaEst{m,it}.hueValFullscreen.delGammaFold{1},gammaEst{m,it}.hueValFullscreen.delGammaFold{2}),4);
                    useXs   = 1:length(useData.val); useXinds= 1:length(useData.val);
                    useXtLabs= useData.val;
                end
                hcols1  = hcols(intersect(0:length(hcols)-1,useData.hueDeg,'stable')+1,:);
                pltData = squeeze(mean(useData.delGamma,1));  % mean across elecs
                pltEst  = squeeze(mean(useEst,1)); 
            end 
            % electrode wise correlation
            for el = 1:size(useEst,1)
                [elCorr{m,pp,it}(el),elCorrSig{m,pp,it}(el)] = corr(squeeze(useData.delGamma(el,:))',squeeze(useEst(el,:))');
            end
        end 
        elCorrM{m,pp} = mean(cell2mat(elCorr(m,pp,:)),3);

        hold(ePlot(m,pp),'on');
        if pp==1
            scatter(ePlot(m,pp),useXs,pltData, [],hcols1);
            plot(ePlot(m,pp),useXs,pltEst,'LineWidth',2,'Color',[0.5 0.5 0.5]);
        else  
            hD = plot(ePlot(m,pp),useXs,squeeze(pltData),'marker','o','linestyle','none');
            hE = plot(ePlot(m,pp),useXs,pltEst,'LineWidth',2,'Color',[0.5 0.5 0.5]);
            set(hD,{'color'},num2cell(hcols1,2)); set(hE,{'color'},num2cell(hcols1,2));
        end 
        text(0.05,0.95,[' n=',num2str(size(useEst,1))],'Color',[0.4 0.4 0.4],'FontWeight','bold','FontSize',fontSizeTicks,'units','normalized','Parent',ePlot(m,pp));
        set(ePlot(m,pp),'XTick',useXs(useXinds),'XTickLabels',round(useXtLabs(useXinds),2),'XTickLabelRotation',45,'ticklength',ticklength,'tickdir','in','fontSize',fontSizeTicks,'fontWeight','bold');    
        xlim(ePlot(m,pp),[useXs(1)-0.5 useXs(end)+0.5]);
        ym=round(max([ym, max(pltData(:))]));
        if m==2, xlabel(ePlot(m,pp),xlabs{pp},'FontWeight','bold','FontSize',fontSizeLarge);
        else,    title(ePlot(m,pp),labs2{pp},'FontWeight','bold','FontSize',fontSizeLarge);
        end
        
        hold(rPlot(m),'on');
        text(0.82,1.20-pp*dd,['n=',num2str(length(elCorr{m,pp}))],'FontWeight','bold','FontSize',fontSizeTicks,'Color',col(pp,:),'units','normalized','Parent',rPlot(m));
        scatter(rPlot(m),repmat(pp,size(elCorrM{m,pp})),elCorrM{m,pp},26,[0.4 0.4 0.4],'filled','Marker','o','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',1);
        mnCor = mean(elCorrM{m,pp}); semCor = std(elCorrM{m,pp})./sqrt(length(elCorrM{m,pp}));
        bar(rPlot(m),pp,mnCor,'FaceColor',col(pp,:),'FaceAlpha',0,'EdgeColor',col(pp,:),'LineWidth',2);
        errorbar(rPlot(m),pp,mnCor,semCor,'Color','k');
    end 
    xlim(rPlot(m),[0.2 5]);  ylim(rPlot(m),[0 1.25]); %    
    set(rPlot(m),'XTick',1:4,'XTickLabel',labs2,'XTickLabelRotation',45,'YTick',0:0.5:1,'YTickLabel',0:0.5:1,'tickdir','in','ticklength',ticklength,'box','off','FontSize',fontSizeTicks,'FontWeight','bold');
    text(-0.34,0.25,'Correlation','Rotation',90,'FontSize',fontSizeSmall,'FontWeight','bold','Units','normalized','Parent',rPlot(m));
    text(-0.46,0.4,['M',num2str(m),'R'],'Rotation',90,'FontSize',fontSizeSmall,'FontWeight','bold','Units','normalized','Parent',rPlot(m));
    text(-0.32,0.15,'Gamma ratio','Rotation',90,'FontSize',fontSizeSmall,'FontWeight','bold','Units','normalized','Parent',ePlot(m,1));
    text(-0.46,0.4,['M',num2str(m),'R'],'Rotation',90,'FontSize',fontSizeSmall,'FontWeight','bold','Units','normalized','Parent',ePlot(m,1));
    for pp=1:4
        ylim(ePlot(m,pp),[0 ym]);
        set(ePlot(m,pp),'YTick',0:2:ym,'YTickLabels',0:2:ym,'YTickLabelRotation',0,'ticklength',ticklength,'tickdir','in');    
    end 
end
text(-0.51,1.1,'A','Rotation',0,'FontSize',fontSizeLabel,'FontWeight','bold','Units','normalized','Parent',ePlot(1,1));
text(-0.40,1.1,'B','Rotation',0,'FontSize',fontSizeLabel,'FontWeight','bold','Units','normalized','Parent',rPlot(1));
set(rPlot(1),'XTickLabel',[]);
if ~showFig
    close(figM);
    close(figH);
end
end

