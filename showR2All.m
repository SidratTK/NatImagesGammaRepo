% showR2 All 
% function for checking independence of features
% this function takes all electrodes of each protocol type and for
% each electrode checks the independence between the variables.
% makes marginal sum, marginal product and svd product based models
% makes figures for typical electrodes.
% plots the distribution across elecs of R2 (variance explained) by model (marginal and Svd based linear model)
% prints out value of separability index of variables for all protocol types
% STK 2023

function [mrgPrd] = showR2All()

folderSourceString=pwd; % where 'data' folder with extracted data is kept

[mrgPrd{1},svdPrd{1},mrgSum{1},subjectNames{1},protocolTypes{1},~,electrodeList{1}] = checkIndependentTuningEnv(folderSourceString,0);
[mrgPrd{2},svdPrd{2},mrgSum{2},subjectNames{2},protocolTypes{2},~,electrodeList{2}] = checkIndependentTuningHueEnv(folderSourceString,0);

titlesP = {{'SF-Ori','Size-Ori','Con-Ori'},{'Hue-Sat','Hue-Val','Sat-Val'}};%cellfun(@(x) x(1:min(6,length(x))),protocolTypes,'un',0); 
tag     = {'','R'};  % 
xTls = {'MSum','MPrd','SVPrd'}; 

figB=figure; colormap jet;
figB.Name = 'Bar Plots of %Variance explained';
figB.Units      = 'centimeters';    figB.PaperType  = 'a4';
figB.PaperUnits = 'centimeters';    figB.PaperSize  = [17 12]; %
figB.Color      = [1 1 1];          figB.PaperOrientation= 'Portrait';
figB.Position   = [0 0 figB.PaperSize]; figB.PaperPosition = [0 0 figB.PaperSize];

plotsR{1}{1} = getPlotHandles(1,2,[0.07 0.65 0.24 0.27],0.03,0.01);
plotsR{1}{2} = getPlotHandles(1,2,[0.41 0.65 0.24 0.27],0.03,0.01);
plotsR{1}{3} = getPlotHandles(1,2,[0.75 0.65 0.24 0.27],0.03,0.01);
plotsR{2}{1} = getPlotHandles(1,2,[0.07 0.15 0.24 0.27],0.03,0.01);
plotsR{2}{2} = getPlotHandles(1,2,[0.41 0.15 0.24 0.27],0.03,0.01);
plotsR{2}{3} = getPlotHandles(1,2,[0.75 0.15 0.24 0.27],0.03,0.01);

fontSizeTicks = 9; fontSizeSmall = 10; fontSizeLarge = 11;
colMs = {[0.2 0.2 0.2],[0.4 0.4 0.4]}; 
col1s = {'b','m','g'}; 
col2s = {[0.4 0.4 0.4],[1 0 0]}; 

disp('Squared correlation R2 represents how much variance is shared between 2 distributions'); 
disp('Or how much variance of A(observed-data) can be explained by B(separable-tuning-estimate)')
disp('Separability index is calculated from SVD decomposition and tells us how well the data can be factored into independent components');

for kk=1:2
    for sub = 1:length(subjectNames{kk}) % 
        numProtocols = length(protocolTypes{kk}); 
        if kk==2 && strcmp(subjectNames{kk}{sub},'kesari'),numProtocols=2; end    % alpa has satval, kesari doesnt
        for pt=1:numProtocols
            hold(plotsR{kk}{pt}(sub),'on');
            ra = cell2mat(cellfun(@(x) x.r2, mrgSum{kk}{sub,pt},'un',0));
            rp = cell2mat(cellfun(@(x) x.r2, mrgPrd{kk}{sub,pt},'un',0));
            rs = cell2mat(cellfun(@(x) x.r2, svdPrd{kk}{sub,pt},'un',0));
            ns=length(rp); 
            meanrs = [mean(ra),mean(rp),mean(rs)]; 
            semrs = [std(ra),std(rp),std(rs)]./sqrt(ns); 
            % pick out Ecog elecs
            if kk==1, ecogInds = electrodeList{kk}{sub,pt}>81;
            else,     ecogInds = false(1,length(electrodeList{kk}{sub,pt}));
            end
            scatter(plotsR{kk}{pt}(sub),  ones(1,sum(~ecogInds)),ra(~ecogInds),20,'o','MarkerFaceColor',col2s{1},'MarkerEdgeColor',col2s{1},'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
            scatter(plotsR{kk}{pt}(sub),2*ones(1,sum(~ecogInds)),rp(~ecogInds),20,'o','MarkerFaceColor',col2s{1},'MarkerEdgeColor',col2s{1},'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
            scatter(plotsR{kk}{pt}(sub),3*ones(1,sum(~ecogInds)),rs(~ecogInds),20,'o','MarkerFaceColor',col2s{1},'MarkerEdgeColor',col2s{1},'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
            scatter(plotsR{kk}{pt}(sub),  ones(1,sum(ecogInds)),ra(ecogInds),20,'o','MarkerFaceColor',col2s{2},'MarkerEdgeColor',col2s{1},'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);
            scatter(plotsR{kk}{pt}(sub),2*ones(1,sum(ecogInds)),rp(ecogInds),20,'o','MarkerFaceColor',col2s{2},'MarkerEdgeColor',col2s{1},'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);
            scatter(plotsR{kk}{pt}(sub),3*ones(1,sum(ecogInds)),rs(ecogInds),20,'o','MarkerFaceColor',col2s{2},'MarkerEdgeColor',col2s{1},'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);
          
            bar(plotsR{kk}{pt}(sub),1,meanrs(1),'FaceAlpha',0,'EdgeColor',col1s{1},'LineWidth',2);
            bar(plotsR{kk}{pt}(sub),2,meanrs(2),'FaceAlpha',0,'EdgeColor',col1s{2},'LineWidth',2);
            bar(plotsR{kk}{pt}(sub),3,meanrs(3),'FaceAlpha',0,'EdgeColor',col1s{3},'LineWidth',2);
            
            errorbar(plotsR{kk}{pt}(sub),1:3,meanrs,semrs,'linestyle','none','color','k','linewidth',1);
            set(plotsR{kk}{pt}(sub),'XTick',1:3,'XTickLabel',xTls,'XTickLabelRotation',90,'YTick',0:0.5:1,'YTickLabel',[],'tickdir','out','box','off','FontSize',fontSizeTicks,'FontWeight','bold');
            xlim(plotsR{kk}{pt}(sub),[0.5 3.5]); ylim(plotsR{kk}{pt}(sub),[0 1.15]); 
            text(0.45,1,['n=',num2str(ns)],'Color','k','FontSize',fontSizeTicks,'FontWeight','bold','Parent',plotsR{kk}{pt}(sub),'units','normalized','color',colMs{sub});
            text(0.02,1,['M',num2str(sub),tag{kk}],'Color','k','FontSize',fontSizeTicks,'FontWeight','bold','Parent',plotsR{kk}{pt}(sub),'units','normalized','color',colMs{sub});
    %        for some significance:
            dosigtest(ra,rp,plotsR{kk}{pt}(sub),[1 2 3], [1 2]); % t-test for whether R2 of MrgPrd is significantly diff from MrgSum
            dosigtest(rp,rs,plotsR{kk}{pt}(sub),[1 2 3], [2 3]); % t-test for whether R2 of MrgPrd is significantly diff from SvdPrd
            if sub==1 
                set(plotsR{kk}{pt}(sub),'YTick',0:0.5:1,'YTickLabel',0:0.5:1);
                text(0.7,1.2,titlesP{kk}{pt},'FontSize',fontSizeLarge,'fontWeight','bold','Parent',plotsR{kk}{pt}(sub),'units','normalized');
                text(-0.45,0.5,sprintf('R^{2}'),'Rotation',90,'fontWeight','bold','fontSize',fontSizeSmall,'Parent',plotsR{kk}{pt}(sub),'units','normalized');
            end
        end
    end 

end
set(plotsR{2}{3}(2),'Visible','off');
end

function dosigtest(x1,x2,hPlot,posX, inds)

if     diff(inds)==1,       ylevel= 1.05;     % where to put significance symbols?
elseif isequal(inds,[1,3]), ylevel= 1.1;
end
xvals = posX(inds) + [0.05 -0.05];     % end points a little inwards
xval1 = [mean(posX(inds))]+[-0.08 0.08];

% [h,p,~,tstats] = ttest(x1,x2);       % paired t test. default sig level =5%
[h,p] = ttest2(x1,x2,'Vartype','unequal');        % unpaired 2sample ttest / Welch test. default sig level =5%
if h==1 && p<=0.005                               % if significantly different distributions, p<0.005 is a subset of p<0.05. 
    plot(hPlot,xvals, [ylevel ylevel],'Color','k');
    plot(hPlot,xval1, [ylevel ylevel],'Color','k','LineStyle','none','Marker','*','MarkerSize',4);
elseif h==1 && p<=0.05                               % if significantly different distributions
    plot(hPlot,xvals, [ylevel ylevel],'Color','k');
    plot(hPlot,mean(posX(inds)), ylevel, 'Color','k','LineStyle','none','Marker','*','MarkerSize',4);
else
    plot(hPlot,xvals, [ylevel ylevel],'Color','k');
    text(mean(posX(inds))-0.2,ylevel-0.03,'ns','fontSize',8,'Parent',hPlot);
end

end

%%%%%%
function [MrgPrd,SVDPrd,MrgSum,subjectNames,protocolTypes,delPowerOut,electrodesOut,varsV,varsH,labsV,labsH] = checkIndependentTuningEnv(folderSourceString,figFlg)
if ~exist('figFlg','var'), figFlg = false; end
dataDir  = 'savedData';
gridType='Microelectrode';

% power related
tBL = [-0.25 0];
tST = [0.25 0.5];
numTapers = 1;
fBandGamma = [30 80];
fBad  = [];

subjectNames = {'alpaH','kesariH'};  
protocolTypes= {'SFOri', 'SizeOri', 'ConOri'};

delPowerOut = cell(length(subjectNames),length(protocolTypes));
electrodesOut=cell(length(subjectNames),length(protocolTypes));
varsV = cell(length(subjectNames),length(protocolTypes));
varsH = cell(length(subjectNames),length(protocolTypes));
labsV = cell(length(subjectNames),length(protocolTypes));
labsH = cell(length(subjectNames),length(protocolTypes));
SVDPrd= cell(length(subjectNames),length(protocolTypes));
MrgPrd= cell(length(subjectNames),length(protocolTypes));
MrgSum= cell(length(subjectNames),length(protocolTypes));
sepIndex=cell(length(subjectNames),length(protocolTypes));
septestP=cell(length(subjectNames),length(protocolTypes));
corrMargSv=cell(length(subjectNames),length(protocolTypes));

for sub = 1:2
    subjectName = subjectNames{sub};

    for pt = 1:length(protocolTypes)
        [expDates, protocolNames, ~, electrodeList, ~] = getProtocolInfoGratings(subjectName,gridType,folderSourceString,protocolTypes{pt},dataDir);
        electrodes = electrodeList{1}; % unwrap from cell
        % get response
        clear delPwrBand meanPwrST meanPwrBL
        if ~strcmp(protocolTypes{pt},'SizeOri') % for full screen
            [delPwrBand,meanPwrST,meanPwrBL,~,~,parameterCombinations,fVals,~,~,~,~,~]= getEnergy(subjectName,expDates{1}{1},protocolNames{1}{1},fBandGamma,fBad,tBL,tST,electrodeList{1},numTapers);
            % ensure these are full screen protocols and not flickering
            sInd = parameterCombinations.sValsUnique==11.5; % fullscreen
            tInd = parameterCombinations.tValsUnique==0;
            parameterCombinations.sValsUnique = parameterCombinations.sValsUnique(sInd);
            parameterCombinations.tValsUnique = parameterCombinations.tValsUnique(tInd);
            delPwrBand = delPwrBand(:,sInd,:,:,:,tInd);
            meanPwrST  = meanPwrST(:,sInd,:,:,:,tInd,:); meanPwrBL = meanPwrBL(:,sInd,:,:,:,tInd,:);
        else
            for e=1:length(electrodes)
                [temp1,temp2,temp3,~,~,pC,fVals,~,~,~,~,~]= getEnergy(subjectName,expDates{1}{e},protocolNames{1}{e},fBandGamma,fBad,tBL,tST,electrodes(e),numTapers);
                sInd = pC.sValsUnique>0.09;   % some alpaH protocol have an extra size. take the ones above 0.09
                pC.sValsUnique = pC.sValsUnique(sInd);
                delPwrBand(e,:,:,:,:,:,:)=temp1(1,sInd,:,:,:,:,:); parameterCombinations = pC;
                meanPwrST(e,:,:,:,:,:,:) =temp2(1,sInd,:,:,:,:,:); meanPwrBL(e,:,:,:,:,:,:)=temp3(1,sInd,:,:,:,:,:);
            end  
        end
        if   strcmp(protocolTypes{pt},'SizeOri'), varV= parameterCombinations.sValsUnique; labV= 'Size';
        elseif strcmp(protocolTypes{pt},'SFOri'), varV= parameterCombinations.fValsUnique; labV= 'SF';
        elseif strcmp(protocolTypes{pt},'ConOri'),varV= parameterCombinations.cValsUnique; labV= 'Cont';
        end
        varV = round(varV,1);
        varH = parameterCombinations.oValsUnique; labH= 'Orientation';
        
        delPowerOut{sub,pt} = delPwrBand;
        electrodesOut{sub,pt} = electrodes;
        varsV{sub,pt} = varV; varsH{sub,pt} = varH;
        labsV{sub,pt} = labV; labsH{sub,pt} = labH;
        
        for e = 1:length(electrodes)
            elecis = electrodes(e);
            delGuse= squeeze(delPwrBand(e,:,:,:,:,:,:));
            stpwr = squeeze(meanPwrST(e,:,:,:,:,:,:));
            blpwr = squeeze(meanPwrBL(e,:,:,:,:,:,:));
            logDelPwr= log10(stpwr./blpwr); 

            if strcmp(subjectName,'alpaH') && elecis==7 && pt==1     % elecis=41 display for typical elec
                   showflg = 1;
            else,  showflg = 0;
            end
            
            [SVDPrd{sub,pt}{e},MrgPrd{sub,pt}{e},MrgSum{sub,pt}{e},sepind,p1,~,margSvdCorr] ...
                = checkIndependentTuning(delGuse,showflg,varH,varV,labH,labV,logDelPwr,fVals,fBandGamma);
            sepIndex{sub,pt}(e) = sepind;     
            septestP{sub,pt}(e) = p1;
            corrMargSv{sub,pt}(e,:) = margSvdCorr;
        end  
    end
end
plotsR2 = showCalculatedResults(protocolTypes,subjectNames,sepIndex,septestP,SVDPrd,MrgPrd,MrgSum,corrMargSv,figFlg);

if figFlg
    fontSizeLabel = 13; fontSizeLarge = 11;
    text(-0.35,0.2,'M 1','Rotation',90,'Parent',plotsR2(1,1),'Units','normalized','FontSize',fontSizeLabel,'FontWeight','bold');
    text(-0.35,0.2,'M 2','Rotation',90,'Parent',plotsR2(2,1),'Units','normalized','FontSize',fontSizeLabel,'FontWeight','bold');
    title(plotsR2(1,1),'SFOri'); title(plotsR2(1,2),'SizeOri'); title(plotsR2(1,3),'ConOri');
    text(-0.25, 0.0,'No. of Electrodes','Rotation',90,'fontWeight','bold','fontSize',fontSizeLarge,'Parent',plotsR2(1,1),'units','normalized');
    text(-0.25, 0.0,'No. of Electrodes','Rotation',90,'fontWeight','bold','fontSize',fontSizeLarge,'Parent',plotsR2(2,1),'units','normalized');
    text(-0.2,-0.2,'Fraction of variance explained','Parent',plotsR2(2,2),'Units','normalized','FontSize',fontSizeLarge,'FontWeight','bold');
end
end
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
    expDates{1}      = {'210817' };  % SFOri FullScreen
    protocolNames{1} = {'GRF_002'};
    expDates{3}      = {'290817' }; % conTFOri
    protocolNames{3} = {'GRF_001'};
    electrodes{2}    =  [82      85      86      88     89];      % ecogs for size protocol
    elecType{2}      = {'ECoG', 'ECoG', 'ECoG', 'ECoG', 'ECoG'};
    expDates{2}      = {'130917' ,'080917' ,'080917' ,'170817' ,'260817'};    % size ori electrode wise
    protocolNames{2} = {'GRF_001','GRF_001','GRF_002','GRF_005','GRF_002'};
    
elseif strcmp(monkeyName,'kesariH')
    expDates{1}      = {'270218' };
    protocolNames{1} = {'GRF_001'};
    expDates{3}      = {'240118' };  % conOri
    protocolNames{3} = {'GRF_002'};
    electrodes{2}    = [85,        86,       88,       89,  ];    % ecogs for size protocol.
    elecType{2}      = {'ECoG',   'ECoG',   'ECoG',   'ECoG' };
    expDates{2}      = {'271217', '030118', '060118', '070118'};   % size ori electrode wise
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
%%%%%
function [MrgPrd,SVDPrd,MrgSum,subjectNames,protocolTypes,delPowerOut,electrodesOut,varsV,varsH,labsV,labsH] = checkIndependentTuningHueEnv(folderSourceString,figFlg)
if ~exist('figFlg','var'), figFlg = false; end
dataDir  = 'savedData';
% power related
tBL = [-0.25 0];
tST = [0.25 0.5];
numTapers = 1;
fBandGamma = [30 80];
fBad  = [];% [48:52, 98:102, 148:152]; % 

% protocol types
protocolTypes= {'HueSatFullscreen','HueValFullscreen','SatVal'};
subjectNames = {'alpa','kesari'};  % hue protocols done for alpa, kesari, not hybrid ecog arrays
numProtocols = [3,2];              % alpa has satval, kesari doesn't

delPowerOut = cell(length(subjectNames),length(protocolTypes));
electrodesOut=cell(length(subjectNames),length(protocolTypes));
varsV = cell(length(subjectNames),length(protocolTypes));
varsH = cell(length(subjectNames),length(protocolTypes));
labsV = cell(length(subjectNames),length(protocolTypes));
labsH = cell(length(subjectNames),length(protocolTypes));
SVDPrd= cell(length(subjectNames),length(protocolTypes));
MrgPrd= cell(length(subjectNames),length(protocolTypes));
MrgSum= cell(length(subjectNames),length(protocolTypes));
sepIndex=cell(length(subjectNames),length(protocolTypes));
septestP=cell(length(subjectNames),length(protocolTypes));
corrMargSv=cell(length(subjectNames),length(protocolTypes));

for sub = 1:length(subjectNames)
        subjectName = subjectNames{sub};
        for pt = 1:numProtocols(sub)    % for pt protocol types
            [expDates, protocolNames, ~, electrodeList, ~] = getProtocolInfoHues(subjectName,protocolTypes{pt},dataDir);
            electrodes = electrodeList{1};
            % get response
            [delPwrBand,meanPwrST,meanPwrBL,~,~,parameterCombinations,fVals,~,~,~,~,~]= getEnergy(subjectName,expDates{1}{1},protocolNames{1}{1},fBandGamma,fBad,tBL,tST,electrodeList{1},numTapers);

            varH = parameterCombinations.oValsUnique; labH = 'Hue';
            if   strcmp(protocolTypes{pt},'HueSize'),            varV = parameterCombinations.sValsUnique; labV = 'Size';
            elseif strcmp(protocolTypes{pt},'HueSatFullscreen'), varV = parameterCombinations.fValsUnique; labV = 'Sat';
            elseif strcmp(protocolTypes{pt},'HueValFullscreen'), varV = parameterCombinations.cValsUnique; labV = 'Val';
            elseif strcmp(protocolTypes{pt},'SatVal'),           varH = parameterCombinations.cValsUnique; labH = 'Val';
                                                                 varV = parameterCombinations.fValsUnique; labV = 'Sat';
            end
            delPowerOut{sub,pt} = delPwrBand;
            electrodesOut{sub,pt} = electrodes;
            varsV{sub,pt} = varV; varsH{sub,pt} = varH;
            labsV{sub,pt} = labV; labsH{sub,pt} = labH;

           % get gamma band response delGuse for each electrode
            for e = 1:length(electrodes)
                elecis = electrodes(e);
                delGuse = squeeze(delPwrBand(e,:,:,:,:,:,:));
                stpwr = squeeze(meanPwrST(e,:,:,:,:,:,:));
                blpwr = squeeze(meanPwrBL(e,:,:,:,:,:,:));
                logDelPwr= log10(stpwr./blpwr);
                
                if (strcmp(subjectName,'kesari') && elecis==74 && pt==1) % strcmp(subjectName,'alpa') && elecis==70  %  display for typical elec
                      showflg = 1;
                else, showflg = 0;
                end

                [SVDPrd{sub,pt}{e},MrgPrd{sub,pt}{e},MrgSum{sub,pt}{e},sepind,p1,~,mrgsvdcorr] ...
                                  = checkIndependentTuning(delGuse,showflg,varH,varV,labH,labV,logDelPwr,fVals,fBandGamma);
                sepIndex{sub,pt}(e) = sepind; 
                septestP{sub,pt}(e) = p1;
                corrMargSv{sub,pt}(e,:) = mrgsvdcorr;
            end 
        end    
end  

plotsR2 = showCalculatedResults(protocolTypes,subjectNames,sepIndex,septestP,SVDPrd,MrgPrd,MrgSum,corrMargSv,figFlg);
if figFlg
    fontSizeLabel = 13; fontSizeLarge = 11;
    text(-0.35,0.2,['M 1','A'],'Rotation',90,'Parent',plotsR2(1,1),'Units','normalized','FontSize',fontSizeLabel,'FontWeight','bold');
    text(-0.35,0.2,['M 2','A'],'Rotation',90,'Parent',plotsR2(2,1),'Units','normalized','FontSize',fontSizeLabel,'FontWeight','bold');
    title(plotsR2(1,1),'Hue Saturation');
    title(plotsR2(1,2),'Hue Value');
    title(plotsR2(1,3),'Saturation Value');
    text(-0.25, 0.0,'No. of Electrodes','Rotation',90,'fontWeight','bold','fontSize',fontSizeLarge,'Parent',plotsR2(1,1),'units','normalized');
    text(-0.25, 0.0,'No. of Electrodes','Rotation',90,'fontWeight','bold','fontSize',fontSizeLarge,'Parent',plotsR2(2,1),'units','normalized');
    text(-0.2,-0.1,'Fraction of variance explained','Parent',plotsR2(2,2),'Units','normalized','FontSize',fontSizeLarge,'FontWeight','bold');
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
    protocolNames{2} = {'GRF_004'}; 
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
    protocolNames{2} = {'GRF_006'}; 
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
        useElecs =   cat(2, useElecs, electrodeList(j));
        useElecInds= cat(2, useElecInds, j);
    end
end          
end

%%%%%
function [SVDPrdOut,MargPrdOut,MargSumOut,sepIndex,pSepInd,hSepInd,margSvCorr] = ...
          checkIndependentTuning(data2DIn,showfig,labH,labV,labH2,labV2,logDelPwr,fVals,fBandGamma)
% this function is to check independence
% data2DIn: input 2D array- Responses (e.g. LFP power in gamma band) to 2 variables like SFOri data.
% computes marginal responses by summing across rows or cols, and SVD singular vectors
% models estimated joint matrix if the 2 distributions were independent- 
% SVDPrdOut - from SVD factors, or as 
% MargPrdOut - product of marginals, or 
% MargSumOut - sum of margs
% sepindex - also calculate separability index using singular vals as in Mazer et al.
% see  Mazer et al., 2002 PNAS who follow a similar analysis.
% the observed joint tuning (data) and estimated separable tuning (prod of
% marginals) are compared via the Pearson correlation coefficient squared,
% which represents the fraction of variance in the data shared by the estimate.
% Note that the product of the marginal tuning may not be the same order of
% magnitude as the observed data, but the correl coeff is not affected by intercept. 
% see Pena & Konishi compare additive (marginals) & multiplicative (factors) models
% use Model fitting to get the separable estimates.
% REFS: Mazer et al., 2002, PNAS; Pena & Konishi, 2001, Science,

% make dummy data example
% data2DIn = rand(4,3);        % 2 dimensional random data OR
% parV = [1:8]; mparV = 4; sparV = 2; margV1 = pdf('Normal',parV,mparV,sparV); margV1 = margV1(:); 
% parH = [1:4]; mparH = 2; sparH = 1; margH1 = pdf('Normal',parH,mparH,sparH);
% data2DIn = margV1 * margH1;
% checkIndependentTuning(data2DIn,1)

if ~exist('showfig','var'), showfig=0;   end

varH   = 1:size(data2DIn,2);
varV   = 1:size(data2DIn,1);
if ~exist('labH','var'), labH = varH; end
if ~exist('labV','var'), labV = varV; end
if ~exist('labH2','var'), labH2 = 'Feature1'; end
if ~exist('labV2','var'), labV2 = 'Feature2'; end

% get marginal sums of row & col
sumVa = repmat(sum(data2DIn,2),[1,size(data2DIn,2)]);
sumHa = repmat(sum(data2DIn,1),[size(data2DIn,1),1]);

% data values tell how strong response a chosen input combination generates. 
% if want the data to behave like probabilities, sum should be 1, so normalise.
% the stimulus combs shown independently, => 2D data is independently sampled. 
% the 3 conditions of being a pdf are met: +ve, sum=1, independent samples 
sumis  = sum(sum(data2DIn));
data2D = data2DIn/sumis;  % subsequent linear corrs hold even if this normalisation not done.

% 1. get marginals (normalised)
margV = sum(data2D,2);    % marginal of Vertical param
margH = sum(data2D,1);    % marginal of Horizontal param

% 2. get SVD factors,
[U,S,V] = svd(data2DIn,'econ');  
S = diag(S); factH = V(:,1); factV = U(:,1);    % val & first singular vectors
sepIndex = S(1)^2/(sum(S.^2));
% 2b. to get significance of sepindex. get random data. 20 times. 
%  ttest tests H0: data in x - y are from 0mean Normal distrib
sibt = zeros(1,20);
for ii=1:20
    a=normrnd(mean(data2DIn(:)),std(data2DIn(:)),size(data2DIn));
    [~,S1,~] = svd(a,'econ');   S1 = diag(S1); 
    sibt(ii) = S1(1)^2/(sum(S1.^2));
end
[hSepInd,pSepInd] = ttest(repmat(sepIndex,size(sibt)),sibt,'alpha',0.005);

% 3. correlation between the factors
[rh,ph] = corr(margH', factH);
[rv,pv] = corr(margV , factV);

% 4A. use model fits for the above marginal model & svd model
% using data2DIn vs data2D or sumH vs margH doesnt change the signficance but only model coeff magnitudes
y    = data2DIn(:);
sumH = sumHa(:);   
sumV = sumVa(:);
d    = dataset(y, sumH, sumV);
fH   = sqrt(S(1))* repmat(factH(:)',[size(data2DIn,1),1]); fH = fH(:);
fV   = sqrt(S(1))* repmat(factV(:) ,[1,size(data2DIn,2)]); fV = fV(:);
d2   = dataset(y, fH, fV);

LMsumMarg = LinearModel.fit(d,'y ~ sumH + sumV - 1'); % -1 to keep dof same.
LMprdMarg = LinearModel.fit(d,'y ~ sumH : sumV');
LMprdSvd  = LinearModel.fit(d2,'y ~ fH : fV');

% 4B. pick up some features from these models here:
r2Svdprd = LMprdSvd.Rsquared.Ordinary;   % R2 reps ratio of variance explained by the model
r2Margprd= LMprdMarg.Rsquared.Ordinary;
r2Margsum= LMsumMarg.Rsquared.Ordinary;
[p_svdprd, f_svdprd, d_svdprd]  = coefTest(LMprdSvd); % p<palpha - model significantly diff than intercept-only model
[p_margprd,f_margprd,d_margprd] = coefTest(LMprdMarg); % f = f value from the f test
[p_margsum,f_margsum,d_margsum] = coefTest(LMsumMarg);

% 5. Update the estimates using the models
data2DIndp = predict(LMprdMarg,[sumH,sumV]); data2DIndp = reshape(data2DIndp,size(data2DIn));
data2DInds = predict(LMsumMarg,[sumH,sumV]); data2DInds = reshape(data2DInds,size(data2DIn));
data2Dsvd  = predict(LMprdSvd ,[fH,fV]);     data2Dsvd  = reshape(data2Dsvd ,size(data2DIn));

% 6. compile outputs
MargPrdOut.Mdl = LMprdMarg; MargPrdOut.r2=r2Margprd; MargPrdOut.fstat=[p_margprd,f_margprd,d_margprd]; 
MargSumOut.Mdl = LMsumMarg; MargSumOut.r2=r2Margsum; MargSumOut.fstat=[p_margsum,f_margsum,d_margsum]; 
SVDPrdOut.Mdl  = LMprdSvd;  SVDPrdOut.r2 =r2Svdprd;  SVDPrdOut.fstat =[p_svdprd, f_svdprd, d_svdprd ]; 
MargPrdOut.Est = data2DIndp; MargPrdOut.Hvec = sumHa(1,:); MargPrdOut.Vvec = sumVa(:,1);
MargSumOut.Est = data2DInds; MargSumOut.Hvec = sumHa(1,:); MargSumOut.Vvec = sumVa(:,1);
SVDPrdOut.Est  = data2Dsvd;  SVDPrdOut.Hvec  = factH';     SVDPrdOut.Vvec = factV;  SVDPrdOut.Svec = S;
margSvCorr     = [rh rv ph pv];

% 7. display
if showfig
    figH = figure; colormap jet;
    figH.Units      = 'centimeters'; figH.PaperType  = 'a4';
    figH.PaperUnits = 'centimeters'; figH.PaperSize  = [11.4 11];% [12.5 10];
    figH.PaperOrientation = 'Portrait'; figH.Color = [1 1 1];
    figH.PaperPosition = [0 0 figH.PaperSize];
    figH.Position = [0 0 figH.PaperSize]; %     figH.PaperUnits = 'normalized';
    
    axp    = getPlotHandles(1,length(varV),[0.14 0.84 0.85 0.11],0.025,0.02); 
    plots1 = getPlotHandles(1,1,[0.79 0.43 0.18 0.22], 0.02,0.02);
    plotMh = getPlotHandles(1,1,[0.23 0.63 0.26 0.10], 0.02,0.02);
    plotMv = getPlotHandles(1,1,[0.06 0.43 0.15 0.18], 0.02,0.02);
    plot2D    = getPlotHandles(1,1,[0.23 0.43 0.26 0.18], 0.02,0.02); 
    plot2Dprd = getPlotHandles(1,1,[0.12 0.08 0.26 0.18], 0.02,0.02); 
    plot2Dsvd = getPlotHandles(1,1,[0.42 0.08 0.26 0.18], 0.02,0.02);
    plot2Dsum = getPlotHandles(1,1,[0.72 0.08 0.26 0.18], 0.02,0.02); 
    caxl = round([min(min(data2D)) max(max(data2D))],2); % data2DIn
    
    hold(plotMh,'on');  hold(plotMv,'on');  axes(plotMv);camroll(90);
    
    dpm = round(max(logDelPwr(:))); dpt = linspace(0,dpm,2);
    for p = 1:length(varV), hold(axp(p),'on');
        line([fBandGamma; fBandGamma],[-2 -2; 2 2],'LineWidth',1,'Color',[0.5 0.5 0.5],'Parent',axp(p));
        plot(axp(p),fVals,zeros(length(fVals)),'LineWidth',1,'Color',[0.5 0.5 0.5]); 
        hp = plot(axp(p),fVals,squeeze(logDelPwr(p,:,:)),'LineWidth',1);
        set(hp,{'color'},num2cell(hsv(length(varH)),2)); 
        xlim(axp(p),[0 130]); ylim(axp(p),[-0.5 dpm]);
        set(axp(p),'XTick',[0 100],'XTickLabel',[0 100],'YTick',dpt,'YTickLabel',[],'YTickLabelRotation',0,'fontSize',9,'fontWeight','bold','tickdir','out','box','off');
        text(0.05,1,[labV2,num2str(labV(p))],'fontSize',9,'fontWeight','bold','Parent',axp(p),'units','normalized');
    end  
   set(axp(1),'YTick',dpt,'YTickLabel',dpt); cols = hsv(length(varH));
   xlabel(axp(3),'Frequency'); 
   text(-0.4,-0.2,'log power','Rotation',90,'Parent',axp(1),'Units','normalized','FontSize',10,'FontWeight','bold'); 
   locs=[2.2 2.4 2.9 3.2 3.7 4 4.6 5];
   for dd = 1:length(varH) 
       text(locs(dd),1.3,num2str(labH(dd)),'color',cols(dd,:),'fontSize',9,'fontWeight','bold','Parent',axp(1),'units','normalized');
   end  

   bar(plots1,1:length(S),100*(sqrt(S.^2)./sum(sqrt(S.^2))),'EdgeColor','none','FaceColor','r');   % magnitude of sing vectors
   xlim(plots1,[0 length(S)]); ylim(plots1,[0 102]);
   text(-0.45,0.5,'%','Rotation',90,'FontSize',10,'FontWeight','bold','Parent',plots1,'Units','normalized');
   text(0.0,-0.30,'S Values','FontSize',10,'FontWeight','bold','Parent',plots1,'Units','normalized');
   set(plots1, 'tickdir','out','YTickLabelRotation',90,'box','off','FontWeight','bold','FontSize',9);
    
   axes(plot2D); imagesc(flip(data2D,1)); clim(plot2D,caxl);
   text(1.20,0.4,labV2,'Rotation',90,'FontSize',10,'FontWeight','bold','Parent',plot2D,'Units','normalized'); % ylabel(plot2D,labV2);
   text(0.15,-0.35,labH2,'FontSize',10,'FontWeight','bold','Parent',plot2D,'Units','normalized');
   c1 = colorbar(plot2D);  c1.Position = [0.57 0.43 c1.Position(3) c1.Position(4)]; c1.Ticks=caxl;
   text(1.5,0.4,' R','Rotation',90,'FontSize',10,'FontWeight','bold','Parent',plot2D,'Units','normalized');
    
   axes(plot2Dprd); imagesc(flip(data2DIndp,1)/sum(sum(data2DIndp))); clim(plot2Dprd,caxl);
   text( 0.15,1.23,'R ~ F*G +1','FontSize',10,'FontWeight','bold','Parent',plot2Dprd,'Units','normalized');
   text(-0.35,0.3,labV2,'Rotation',90,'FontSize',10,'FontWeight','bold','Parent',plot2Dprd,'Units','normalized');
   text(0.3,1.07, ['R2=',num2str(round(r2Margprd,2))],'Parent',plot2Dprd,'Units','normalized','FontSize',9,'FontWeight','bold');
    
   axes(plot2Dsum); imagesc(flip(data2DInds,1)/sum(sum(data2DInds))); clim(plot2Dsum,caxl); % xlabel(plot2Dsum,labH2); % ylabel(plot2Dsum,labV2);
   text( 0.15,1.23,'R ~ F+G','FontSize',10,'FontWeight','bold','Parent',plot2Dsum,'Units','normalized'); % remove 1+, 081121
   text(0.3,1.07, ['R2=',num2str(round(r2Margsum,2))],'Parent',plot2Dsum,'Units','normalized','FontSize',9,'FontWeight','bold');

   axes(plot2Dsvd); imagesc(flip(data2Dsvd,1)/sum(sum(data2Dsvd))); clim(plot2Dsvd,caxl);
   text(0.1,1.23,'R~ U*s1*V +1','FontSize',10,'FontWeight','bold','Parent',plot2Dsvd,'Units','normalized');
   text(0.15,-0.35,labH2,'FontSize',10,'FontWeight','bold','Parent',plot2Dsvd,'Units','normalized');
   text(0.3,1.07, ['R2=',num2str(round(r2Svdprd,2)) ],'Parent',plot2Dsvd,'Units','normalized','FontSize',9,'FontWeight','bold');

   plot(plotMh,varH,margH,'LineWidth',2,'Color','b'); 
   plot(plotMv,varV,margV,'LineWidth',2,'Color','b'); 
   plot(plotMh,varH,factH/sum(factH),'LineWidth',2,'Color','r','Linestyle','--');
   plot(plotMv,varV,factV/sum(factV),'LineWidth',2,'Color','r','Linestyle','--'); %   xlabel(plotMh,labH2);
   text(0.9,0.9,'Marginal G','Rotation',0,'Parent',plotMh,'Units','normalized','FontSize',9,'FontWeight','bold','Color','b');
   text(0.9,0.6,'S.vector V','Rotation',0,'Parent',plotMh,'Units','normalized','FontSize',9,'FontWeight','bold','Color','r');
   text(0.07,0.9,'Marginal F','Rotation',90,'Parent',plotMv,'Units','normalized','FontSize',9,'FontWeight','bold','Color','b');
   text(0.18,0.9,'S.vector U','Rotation',90,'Parent',plotMv,'Units','normalized','FontSize',9,'FontWeight','bold','Color','r');

   indsH = 1:length(varH); if length(varH)>4, indsH = 1:2:length(varH); end
   set(plot2D   ,'XTick',varH(indsH),'XTickLabel',labH(indsH),'YTick',varV,'YTickLabel',flip(labV),'yaxisLocation','right','XTickLabelRotation',0,'tickdir','out','box','on','fontWeight','bold','FontSize',9);
   set(plot2Dprd,'XTick',varH(indsH),'XTickLabel',labH(indsH),'YTick',varV,'YTickLabel',flip(labV),'XTickLabelRotation',0,'tickdir','out','box','off','fontWeight','bold','FontSize',9)
   set(plot2Dsum,'XTick',varH(indsH),'XTickLabel',labH(indsH),'YTick',varV,'YTickLabel',[],'XTickLabelRotation',0,'tickdir','out','box','off','fontWeight','bold','FontSize',9)
   set(plot2Dsvd,'XTick',varH(indsH),'XTickLabel',labH(indsH),'YTick',varV,'YTickLabel',[],'XTickLabelRotation',0,'tickdir','out','box','off','fontWeight','bold','FontSize',9)
   
   set(plotMh,'XTick',varH(indsH),'XTickLabel',[],'YTick',0:0.15:0.3,'YTickLabel',0:0.15:0.3,'tickdir','out','box','off','fontWeight','bold','FontSize',9);  % labH(indsH)
   set(plotMv,'XTick',varV,'XTickLabel',[],'YTick',0:0.15:0.3,'YTickLabel',0:0.15:0.3,'YTickLabelRotation',90,'tickdir','out','box','off','fontWeight','bold','FontSize',9);  
   xlim(plotMh,[varH(1)-0.5 varH(end)+0.5]); xlim(plotMv,[varV(1)-0.5 varV(end)+0.5]); 
   ylim(plotMh,[0 0.35]); ylim(plotMv,[0 0.35]);
   
   text(-0.75,1.3,'A','Parent',axp(1),'Units','normalized','FontSize',12,'FontWeight','bold');
   text(-0.80,1.1,'B','Parent',plotMh,'Units','normalized','FontSize',12,'FontWeight','bold');
   text(-0.5,1.30,'C','Parent',plots1,'Units','normalized','FontSize',12,'FontWeight','bold');
   text(-0.35,1.3,'D','Parent',plot2Dprd,'Units','normalized','FontSize',12,'FontWeight','bold');

    disp(['Separability index from SVD factors is ',num2str(sepIndex)]);
    disp(['Variance explained by svd product model: ',num2str(100*r2Svdprd),'% .',...
        ' model F-val= ',num2str(round(f_svdprd,2)),' p= ',num2str(round(p_svdprd,3))]);
    disp(['Variance explained by marginals product model: ',num2str(100*r2Margprd),'% .',...
        ' model F-val= ',num2str(round(f_margprd,2)),' p= ',num2str(round(p_margprd,3))]);
    disp(['Variance explained by marginals sum model: ',num2str(100*r2Margsum),'% .',...
        ' model F-val= ',num2str(round(f_margsum,2)),' p= ',num2str(round(p_margsum,3))]);
end


end

function [plotsR2] = showCalculatedResults(protocolTypes,subjectNames,sepIndex,septestP,SVDPrd,MrgPrd,MrgSum,corrMargSv,figFlg)
if ~exist('figFlg','var'), figFlg = false; end

palph = 0.05;
plotsR2=[];

if figFlg
    figH=figure;
    figH.Name = 'Histogram of variance explained by model across electrodes';
    figH.Units      = 'centimeters';    figH.PaperType  = 'a4';
    figH.PaperUnits = 'centimeters';    figH.PaperSize  = [17 10]; %
    figH.Color      = [1 1 1];          figH.PaperOrientation= 'Portrait';
    figH.Position   = [0 0 figH.PaperSize]; figH.PaperPosition = [0 0 figH.PaperSize];
    plotsR2 = getPlotHandles(2,3,[0.1 0.1 0.87 0.8],0.06,0.1);
    fontSizeTicks = 9; fontSizeSmall = 10;
end

for sub = 1:2  
    for pt=1:3
        if isempty(SVDPrd{sub,pt}), continue;  
        end
        sepind = round(mean(sepIndex{sub,pt}),2);ss = round(std(sepIndex{sub,pt}),2);
        hofsep = sum(septestP{sub,pt}<=palph);
  
        rSvd = cell2mat(cellfun(@(x) x.r2, SVDPrd{sub,pt},'un',0));
        r2svd = round(mean(rSvd),2); sr2svd = round(std(rSvd),2);
        n     = length(rSvd);
        pfsvd = cell2mat(cellfun(@(x) x.fstat(1), SVDPrd{sub,pt},'un',0));
        pfsvd = sum(pfsvd <= palph);  % p from F test
        
        rPrd = cell2mat(cellfun(@(x) x.r2, MrgPrd{sub,pt},'un',0));
        r2prd = round(mean(rPrd),2); sr2prd = round(std(rPrd),2);
        pfprd = cell2mat(cellfun(@(x) x.fstat(1), MrgPrd{sub,pt},'un',0));
        pfprd = sum(pfprd <= palph);
        
        rSum = cell2mat(cellfun(@(x) x.r2, MrgSum{sub,pt},'un',0));
        r2sum = round(mean(rSum),2); sr2sum = round(std(rSum),2);
        pfsum = cell2mat(cellfun(@(x) x.fstat(1), MrgSum{sub,pt},'un',0));
        pfsum = sum(pfsum <= palph);
         
        disp(' -- ');
        disp(['  ', protocolTypes{pt},'  ',subjectNames{sub},' :',',  n = ',num2str(n)]);
        disp([' Separability index (svd) = ',num2str(sepind),'+-',num2str(ss),'SD',' ttest h=0 for ',num2str(hofsep)]);
        disp([' SVD r2 = ',num2str(r2svd),'+-',num2str(sr2svd), 'SD. Significant F-stat p-val for ',num2str(pfsvd),' of ',num2str(n)]);
        disp(['MrgP r2 = ',num2str(r2prd),'+-',num2str(sr2prd), 'SD. Significant F-stat p-val for ',num2str(pfprd),' of ',num2str(n)]);
        disp(['MrgS r2 = ',num2str(r2sum),'+-',num2str(sr2sum), 'SD. Significant F-stat p-val for ',num2str(pfsum),' of ',num2str(n)]);

        inds = (corrMargSv{sub,pt}(:,3)<=palph) & (corrMargSv{sub,pt}(:,4)<=palph);
        mr = mean(abs([corrMargSv{sub,pt}(:,1);corrMargSv{sub,pt}(:,2)])); % horiz and vertical vecs %corrMargSv{sub,pt}(inds,1)
        sr = std(abs([corrMargSv{sub,pt}(:,1);corrMargSv{sub,pt}(:,2)]));
        disp(['Corr bw SVD factors and marginals was significant for ',num2str(sum(inds)),' of ',num2str(n)]);
        disp(['Avg mag of Corr coef =', num2str(round(mr,2)),' +-',num2str(round(sr,2)),' SD']);
        
        [h,p,~,~] = ttest2(rPrd,rSum,'Vartype','unequal'); % t-test for whether R2 of MrgPrd is significantly diff from MrgSum
        disp(['Unpaired ttest2 on R2 of Marg Sum & Prd. H = ',num2str(h),'P = ',num2str(round(p,3))]);
        [h,p,~,~] = ttest2(rPrd,rSvd,'Vartype','unequal'); % t-test for whether R2 of MrgPrd is significantly diff from svd?
        disp(['Unpaired ttest2 on R2 of Marg Prd & svd. H = ',num2str(h),'P = ',num2str(round(p,3))]);
        
        if figFlg % plot the distribution of R2 of marg product model
            hold(plotsR2(sub,pt),'on');
            hi = histogram(plotsR2(sub,pt),rPrd,0:0.05:1);
            hi.FaceAlpha = 0.6; hi.EdgeColor = 'none'; hi.FaceColor = 'b';
            text(0.05,0.95,['n= ',num2str(length(rPrd))],'Color','b','FontSize',fontSizeSmall,'FontWeight','bold','Parent',plotsR2(sub,pt),'units','normalized');
            set(plotsR2(sub,pt),'XTick',0:0.5:1,'XTickLabel',[0,0.5,1],'YTick',0:5:35,'YTickLabel',0:5:35,'tickdir','out','box','off','FontSize',fontSizeTicks,'FontWeight','bold');
        end 
    end
end

end