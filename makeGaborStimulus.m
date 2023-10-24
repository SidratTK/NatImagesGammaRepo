% makeGaborStimulus(gaborStim,aVals,eVals)
% G(x,y) = exp(-((x'^2) + (y'^2))/2) sin(2*pi*f*x')
% where 
% x' = xcos(theta)-ysin(theta)
% y' = xsin(theta)+ycos(theta)

% Adding option to draw a colored luminance grating specified by gaborStim.hueDeg (0 to 360) and gaborStim.sat(0 to 1).
function [gaborPatch,aperature] = makeGaborStimulus(gaborStim,aVals,eVals,showGabor)

if ~exist('showGabor','var');               showGabor=0;                end

azi = gaborStim.azimuthDeg;
ele = gaborStim.elevationDeg;
sf  = gaborStim.spatialFreqCPD;
ori = gaborStim.orientationDeg;
C   = gaborStim.contrastPC/2;

if isfield(gaborStim,'spatialFreqPhaseDeg')
    phi = gaborStim.spatialFreqPhaseDeg*(pi/180);
else
    phi=0;
end

if length(gaborStim.radiusDeg)==1    % Gabor
    radMax = gaborStim.radiusDeg;
    radMin = 0;
else
    radMax = gaborStim.radiusDeg(2);
    radMin = gaborStim.radiusDeg(1);
end

% Make a 2D grating
theta = pi*ori/180; % converting to radians

nE=length(eVals);
nA=length(aVals);

grating = zeros(nE,nA);
aperature = zeros(nE,nA);

for e=1:nE
    for a=1:nA
        xg = aVals(a)*cos(theta) - eVals(e)*sin(theta);
        grating(e,a) = C*sin(2*pi*sf*xg + phi);
        
        distance = sqrt((aVals(a)-azi)^2+(eVals(e)-ele)^2);
        
        if (distance<=radMax) && (distance>=radMin)
            aperature(e,a) = 1;
        else
            aperature(e,a) = 0;
        end
    end
end

% Gaussian
params(1) = azi;
params(2) = ele;
params(3) = gaborStim.sigmaDeg;
params(4) = params(3);
params(5) = 0;
params(6) = 1;
[~,GaussianEnvelope,boundaryX,boundaryY] = gauss2D(params,aVals,eVals,[]);

% set everything outside radius to zero
gaborPatch = (50+GaussianEnvelope.*aperature.*grating)/100; % Now between 0 and 1 to maintain uniformity with colored gratings. Earlier was between 0 and 100.

% If hue and sat values are provided, fill with color
if isfield(gaborStim,'hueDeg')
    gaborHSV(:,:,3) = gaborPatch; % Val between 0-1
    gaborHSV(:,:,1) = gaborStim.hueDeg/360; % Hue between 0-1
    gaborHSV(:,:,2) = gaborStim.sat; % Sat between 0-1
    
    gaborRGB = hsv2rgb(gaborHSV); % Convert from HSV to RGB
    
    % setting background back to gray 0.5
    gaborRGB(:,:,1) = (~aperature)*0.5+gaborRGB(:,:,1).*(aperature);
    gaborRGB(:,:,2) = (~aperature)*0.5+gaborRGB(:,:,2).*(aperature);
    gaborRGB(:,:,3) = (~aperature)*0.5+gaborRGB(:,:,3).*(aperature);

    gaborPatch=gaborRGB;
end

% Note that previous we were using pcolor for which the y-axis goes from
% small to large values like the monitor. But now we have to use imagesc to
% show color images, for which the y-axis orientation gets flipped.
if showGabor
    colormap('gray');
    subplot(221);
    imagesc(aVals,eVals,grating); colorbar;
    subplot(222);
    imagesc(aVals,eVals,GaussianEnvelope); colorbar;
    subplot(223);
    imagesc(aVals,eVals,aperature); colorbar;
    subplot(224);
    imagesc(aVals,eVals,gaborPatch); colorbar
    hold on;
    plot(boundaryX,boundaryY,'r');
end