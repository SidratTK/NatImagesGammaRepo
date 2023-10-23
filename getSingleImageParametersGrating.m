% function to check for grating like features in an image patch
% inputs: imageHSV. 3d image in HSV format. 
%         only the V (3rd) layer is used. ranges from 0-1
%         imageAxesDeg: x and y axis values in degrees
%         rfCenterDeg: location around which patch is checked
%         useSigmasDeg: sizes. sigmas of patches to be checked
% outputs: gaborParams: the best approximation of the image patch in terms of a Gabor.
%          Has fields spatialFreqCPD, orientationDeg, sigmaDeg, radiusDeg, spatialFreqPhaseDeg
%          Has a field categoryGrating which flags if the patch qualifies
%          for a grating feature. Should have SF within a range (not too
%          low, not too high) and Orientation Variance above a threhsold.
%          also saved Orientation Variance at the chosen SF. 
%          

function [gaborParams] = getSingleImageParametersGrating(imageHSV,imageAxesDeg,rfCenterDeg,useSigmasDeg)
if ~exist('useSigmasDeg','var'), useSigmasDeg = 0.3:0.3:2.1;    end

imgVal = imageHSV(:,:,3);      % val layer
ovThresh = 20;
sfLims   = [0.4 8];

gaborParams = getGaborParams(imgVal,imageAxesDeg,rfCenterDeg,useSigmasDeg);
% put some conditions:
if gaborParams.spatialFreqCPD>=sfLims(1) && gaborParams.spatialFreqCPD<=sfLims(2) ... % if SF is too low or too high, not grating
   && gaborParams.oriVar>= ovThresh
      gaborParams.categoryGabor = true; 
else, gaborParams.categoryGabor = false;
end
end

function gaborStim = getGaborParams(imgVal,imageAxesDeg,rfCenterDeg,useSigmasDeg)
% initialise a Gabor
gaborStim.azimuthDeg  = rfCenterDeg(1);
gaborStim.elevationDeg= rfCenterDeg(2);
gaborStim.contrastPC  = 100;
gaborStim.spatialFreqPhaseDeg=-90; 
gaborStim.sigmaDeg      = 0;
gaborStim.orientationDeg= 0;
gaborStim.spatialFreqCPD= 0;
gaborStim.radiusDeg     = 3*gaborStim.sigmaDeg;
gaborStim.oriVar        = 0;

imageDFT= cell(1,length(useSigmasDeg));
peakProdNorm = cell(1,length(useSigmasDeg)); normOfFilt=cell(1,length(useSigmasDeg));
peakProd= cell(1,length(useSigmasDeg));
peakSfs = cell(1,length(useSigmasDeg));
peakOrs = cell(1,length(useSigmasDeg));
oriVars = cell(1,length(useSigmasDeg));

for rr= 1:length(useSigmasDeg) % each size
    checkSigmaDeg = useSigmasDeg(rr);
    [imageDFT{rr},Fx,Fy] = getSpectraSingleImageFun(imgVal,rfCenterDeg,checkSigmaDeg,imageAxesDeg);
    [peakSfs{rr},peakOrs{rr},peakProd{rr},~]= getImageDftPeaks(imageDFT{rr},Fx,Fy);
    peakOrs{rr}  = peakOrs{rr}(peakSfs{rr}>=0.1); % v.low sf, its jus due to nat img stat & may mask real gratng feature bcz of hi prod
    peakProd{rr} = peakProd{rr}(peakSfs{rr}>=0.1);
    peakSfs{rr}  = peakSfs{rr}(peakSfs{rr}>=0.1);

    if ~isempty(peakSfs{rr})
        for ii=1:length(peakSfs{rr})
            % normalize the product by the self product of the Gabor. Larger sizes have a larger product. 
            gaborStim.sigmaDeg = useSigmasDeg(rr);
            gaborStim.radiusDeg= 3*gaborStim.sigmaDeg;
            gaborStim.orientationDeg= peakOrs{rr}(ii);
            gaborStim.spatialFreqCPD= peakSfs{rr}(ii);
            gaborIs = makeGaborStimulus(gaborStim,imageAxesDeg.xAxisDeg,imageAxesDeg.yAxisDeg,0);
            gaborIs = gaborIs-0.5;   % get from [0 1] to [-0.5 0.5]
            normOfFilt{rr}(ii) = sqrt(sum(sum(gaborIs.*gaborIs))); % euclidean norm of filt   
    
            % calculate the Ori variance at this SF?
            oriVars{rr}(ii) = getOriVarFrom2D(imageDFT{rr},Fx,Fy,peakSfs{rr}(ii),normOfFilt{rr}(ii)); % 

        end
        % which is the peak sf, ori at this Size?
        maxProdInd = find(peakProd{rr}==max(peakProd{rr}));
        if length(maxProdInd)>1 
            maxProdInd = maxProdInd(oriVars{rr}(maxProdInd)==max(oriVars{rr}(maxProdInd)));
        end
        peakProd{rr} = peakProd{rr}(maxProdInd(1));  % use the first one, in case same product. mirror image
        peakSfs{rr}  = peakSfs{rr}(maxProdInd(1));
        peakOrs{rr}  = peakOrs{rr}(maxProdInd(1));
        oriVars{rr}  = oriVars{rr}(maxProdInd(1));
        normOfFilt{rr}=normOfFilt{rr}(maxProdInd(1));
        peakProdNorm{rr}= peakProd{rr}./normOfFilt{rr};
    end
end 
% check for max among sizes
if any(cell2mat(cellfun(@(x) ~isempty(x), peakProdNorm,'un',0)))
    useInd = find(cell2mat(cellfun(@(x) ~isempty(x), peakProdNorm,'un',0)));
    peakProdNormMat = cell2mat(peakProdNorm);
    maxProdInd = peakProdNormMat==max(peakProdNormMat);
    maxProdInd = useInd(maxProdInd);
    gaborStim.sigmaDeg      = useSigmasDeg(maxProdInd(1));
    gaborStim.orientationDeg= peakOrs{maxProdInd(1)};
    gaborStim.spatialFreqCPD= peakSfs{maxProdInd(1)};
    gaborStim.radiusDeg     = 3*gaborStim.sigmaDeg;
    gaborStim.oriVar        = oriVars{maxProdInd(1)};
    
    % check for phase
    gaborStim = getPhase(imgVal,gaborStim,imageAxesDeg);
    % check for contrast 
    gaborStim = getContrast(imgVal,gaborStim,imageAxesDeg); 
    
end
end

function [imageDFT,Fx,Fy] = getSpectraSingleImageFun(imgVal,rfCenterDeg,checkSigmaDeg,imageAxesDeg)
        imgVal = imgVal-(mean(imgVal(:)));   % remove dc if any? 
        % gaussian mask to only keep the features at center & fade the rest with smooth edges.
        params= [rfCenterDeg checkSigmaDeg checkSigmaDeg 0 1 0]; % center azi ele, sd major minor, angle of rotation wrt coord system, scaling factor, additive constant
        [~,gaussianEnvelope,~,~] = gauss2D(params,imageAxesDeg.xAxisDeg,imageAxesDeg.yAxisDeg,[]);
        imgValMasked = imgVal .* gaussianEnvelope;
        [imageDFT,Fx,Fy] = getImageDFT(imgValMasked,imageAxesDeg);  % get 2D fourier xfm of image layers & its 1D radial avg.     
end  
function [imageDFT,Fx,Fy]= getImageDFT(imageIn,imageAxesDeg)
        [hIm,wIm] = size(imageIn); % 2D image
        % get frequency domain axes:
        dFy = 1/(imageAxesDeg.yAxisDeg(end)- imageAxesDeg.yAxisDeg(1));  % inv of spatial spread (deg) 
        dFx = 1/(imageAxesDeg.xAxisDeg(end)- imageAxesDeg.xAxisDeg(1));  % is freq domain resolution (per deg)
        Fs_x= wIm.*dFx;               Fs_y= hIm.*dFy;     % image pixels per degree - sampling freq
        Fx  = -Fs_x/2:dFx:Fs_x/2-dFx; Fy  = -Fs_y/2:dFy:Fs_y/2-dFy;
        imageDFT = fftshift(fft2(imageIn));                % 2D FFT
end 
function [peakSfs,peakOrs,peakPrds,peakIndsPsd]= getImageDftPeaks(imageDFT,Fx,Fy)
% checks the 2d image fourier transform for peaks 
imagePwr = (abs(imageDFT));
imagePwr = imagePwr./max(max(imagePwr)); % normalise by maximum
% find peaks. 
peaksInCols = islocalmax(imagePwr,1,'FlatSelection','center','MinProminence',0.1);
peaksInRows = islocalmax(imagePwr,2,'FlatSelection','center','MinProminence',0.1);
peakIndsPsd = peaksInCols & peaksInRows;
[rowind,colind]=find(peakIndsPsd);
% find the sf and ori corresp to the peaks 
peakSfs=zeros(1,length(rowind));peakOrs=zeros(1,length(rowind));peakPrds=zeros(1,length(rowind));
Fy = flip(Fy); % to match with the y axis of 2d matrix 
for pknum = 1:length(rowind)
    peakSfs(pknum) = sqrt(Fy(rowind(pknum))^2 + Fx(colind(pknum))^2);
    peakOrs(pknum) = atan2d(Fy(rowind(pknum)) , Fx(colind(pknum)));
    peakPrds(pknum)= abs(imageDFT(rowind(pknum),colind(pknum)));
end
end  
function [oriVar] = getOriVarFrom2D(imageDFT,Fx,Fy,sf,normFactor)
if ~exist('normFactor','var'), normFactor=1; end
imageFAmp= abs(imageDFT);
[X, Y] = meshgrid(Fx, flip(Fy));  % flip Fy so -ve is lower quadrants as for imageDFT 
sfsAll = nan(size(imageDFT));
angD   = nan(size(imageDFT));
for rows=1:size(imageDFT,1)
    for cols=1:size(imageDFT,2)
        sfsAll(rows,cols)= sqrt( (X(rows,cols))^2 + (Y(rows,cols))^2 );
        angD(rows,cols)= atan2d( Y(rows,cols),X(rows,cols) );
    end 
end 
useInds1= sfsAll==sf;
useProds= imageFAmp(useInds1);% all ori responses at this SF
oriVar  = var(useProds./normFactor);     
end
function [gaborParams]= getPhase(imageVal,gaborParams,imageAxesDeg)
% initialise gabor inputs for making a gabor 
gaborStim = gaborParams;
gaborStim.contrastPC = 100;
usePhases = [0 90]; %  get prod at 2 phases 
useProd   = zeros(size(usePhases));
for pp=1:length(usePhases)
    gaborStim.spatialFreqPhaseDeg=usePhases(pp);
    gaborIs = makeGaborStimulus(gaborStim,imageAxesDeg.xAxisDeg,imageAxesDeg.yAxisDeg,0);
    useProd(pp) = (sum(sum((gaborIs-0.5).*(imageVal-mean(imageVal(:))))));
end 
phout1 = atan2d(useProd(2),useProd(1));
gaborParams.spatialFreqPhaseDeg = phout1;
end
function [gaborStim]= getContrast(imageVal,gaborParams,imageAxesDeg)
gaborStim = gaborParams;
% what is the Mich contrast in patch?
gaborParams.radiusDeg = gaborParams.sigmaDeg;
[~,ellipsePixels] = makeGaborStimulus(gaborParams,imageAxesDeg.xAxisDeg,imageAxesDeg.yAxisDeg); % generate a circular mask      
[rowinds,colinds] = find(ellipsePixels);
numinds = length(rowinds);
chosenPixVal = zeros(1,numinds);
for ind = 1:numinds
    chosenPixVal(ind) =  imageVal(rowinds(ind),colinds(ind));
end
maxLum= max(chosenPixVal);
minLum= min(chosenPixVal);

conMich = (maxLum-minLum)/(maxLum+minLum);
gaborStim.contrastPC = conMich*100;
end

