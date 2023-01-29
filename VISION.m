% VISION: VIdeo StabilisatION using automatic features selection
% ==============================================================
% 
% VISION.m is a Matlab function aiming at stabilising videos for image 
% velocimetry analyses in rivers. It is a command line function without GUI
% at the moment. Written in Matlab R2020a. Computer Vision Toolbox is required. 
% 
% 
% The function has one input and six configurable settings:
% --------------------------------------------------------
% 
% *************************************************************************
% * INPUT:                                                                *  
% *                                                                       *
% * [1] filename: The name of the video file, e.g.: 'RawVideo.avi'        *
% *                                                                       *
% *************************************************************************
%  
% ********************************************************************************************************
% * CONFIGURABLE SETTINGS:                                                                               *
% *                                                                                                      *
% * [2] NameStabilisedVideo: Name of the stabilised video, e.g.: 'StabilisedVideo'.                      * 
% *                          The stabilised video is saved by default in the avi format.                 * 
% *                          Default: StabilisedVideo + Feature detection algorithms.                    *
% *                                                                                                      *
% * [3] algorithms: The algorithm(s) used for features' detection. At the moment, a maximum              *
% *                 of two algorithms can be chosen given to the user flexibility within the             *
% *                 analysis. The algorithms must be written within parentheses,                         *
% *                 e.g.: {'FAST','SURF'} or {'FAST'}.                                                   *
% *                 The following algorithms have been implemented:                                      *
% *                 'FAST' / 'MINEIGEN' / 'HARRIS' / 'SURF' / 'KAZE' / 'BRISK' / 'ORB' /                 *
% *                 Default: SURF.                                                                       *
% *                                                                                                      *
% * [3] PercentualStrongestFeatures: Percentual value of the strongest features detected to be           *
% *                                  used for stabilisation analyses. Range of accepted values:          *
% *                                  ]0, 100]. 100 means all the features to be considered.              *
% *                                  An uniform filter is afterwards applied to remain with the          *
% *                                  50% of the value entered by the user of uniformly distributed       *
% *                                  features. Default: 100%.                                            *
% *                                                                                                      *
% * [4] TransformType: Transformation type among 'similarity', 'affine', or 'projective'.                *
% *                    Minimum number of matched points: similarity < affine < projective.               *
% *                    Better accuracy of estimated transformation when greater the number               *
% *                    of matched points is. Default: similarity.                                        *
% *                                                                                                      *
% * [5] ROISelection: Two possible options: i) Binary decision, and ii) Rectangular ROI introduced       *
% *                                         a priori in a vector format. Binary decision giving the      *
% *                                         possibilitity to define a ROI for stabilisation analysis.    *
% *                                         1 means ROI, while 0 means all the image in question.        *
% *                                         Rectangular ROI in a vector format: e.g., [0 0 1000 1000].   *
% *                                         Default: all the field of view.                              *
% *                                                                                                      *
% * [6] StabilisationViewer: Binary decision giving the possibilitity to open a viewer to see how the    *
% *                          stabilisation goes. 1 means viewer, while 0 no viewer. Default: 1 (viewer)  *
% *                                                                                                      *
% ********************************************************************************************************
%  
% *************************************************************************
% * OUTPUTS:                                                              *  
% *                                                                       *
% * [1] Stabilised video                                                  *
% * [2] Number of Frames (NFrames)                                        *
% * [3] Frame Per Seconds (FPS)                                           *
% * [4] Region Of Interest (ROI)                                          *
% *                                                                       *
% *************************************************************************
% 
% Examples of how to call VISION:
% ------------------------------
% 
% [NFrames,FPS,ROI] = VISION('Belgrade_15frames.avi')
% [NFrames,FPS,ROI] = VISION('Belgrade_15frames.avi','NameStabilisedVideo','TestVideo')
% [NFrames,FPS,ROI] = VISION('Belgrade_15frames.avi','NameStabilisedVideo','TestVideo','algorithms',{'FAST'})
% [NFrames,FPS,ROI] = VISION('Belgrade_15frames.avi','NameStabilisedVideo','TestVideo','algorithms',{'FAST'},'PercentualStrongestFeatures',30)
% [NFrames,FPS,ROI] = VISION('Belgrade_15frames.avi','NameStabilisedVideo','TestVideo','algorithms',{'FAST'},'PercentualStrongestFeatures',30,'TransformType','affine')
% [NFrames,FPS,ROI] = VISION('Belgrade_15frames.avi','NameStabilisedVideo','TestVideo','algorithms',{'FAST'},'PercentualStrongestFeatures',30,'TransformType','affine','ROISelection',[1 1 1000 1000])
% [NFrames,FPS,ROI] = VISION('Belgrade_15frames.avi','NameStabilisedVideo','TestVideo','algorithms',{'FAST'},'PercentualStrongestFeatures',30,'TransformType','affine','ROISelection',[1 1 1000 1000],'StabilisationViewer',0)
% 
% 
% This script was prepared by Alonso Pizarro (1), Silvano F. Dal Sasso (2) & Salvatore Manfreda (3)
% 
% (1) Universidad Diego Portales, Chile | alonso.pizarro[at]mail.udp.cl (https://orcid.org/0000-0002-7242-6559)
% (2) Università Degli Studi Della Basilicata, Italy | silvano.dalsasso[at]unibas.it (https://orcid.org/0000-0003-1376-7764)
% (3) University of Naples Federico II, Italy | salvatore.manfreda[at]unina.it (https://orcid.org/0000-0002-0225-144X)
% 15/02/2021 - V0.0.1
% 08/04/2021 - V0.0.2
% 01/03/2021 - V0.0.3
% 
% Copyright (C) 2021 Alonso Pizarro, Silvano F. Dal Sasso & Salvatore Manfreda
% This program is free software (BSD 3-Clause) and distributed WITHOUT ANY
% WARRANTY.


function [NFrames,FPS,ROI] = VISION(filename, varargin)
%% Checking inputs and Setting defaults values
if nargin < 1
    error('Not enugh input arguments')
elseif nargin > 11
    error('Too many input arguments')    
end

ip = inputParser;
ip.CaseSensitive = true;   
    
% required input arguments
addRequired(ip, 'filename', @(filename) ischar(filename)) 

% optional input arguments
addParameter(ip, 'NameStabilisedVideo', 'StabilisedVideo', @ischar) %Name of stabilised video. Default: StabilisedVideo + Feature detection algorithms
addParameter(ip, 'algorithms', {'SURF'}, @iscell)                   %Features detection algorithms. Default: SURF
addParameter(ip, 'PercentualStrongestFeatures', 100, @isnumeric)    %Percentual strongest features. Default: 100%
addParameter(ip, 'TransformType', 'similarity', @ischar)            %Transform type. Default: similarity
addParameter(ip, 'ROISelection', 0, @isnumeric)                     %ROI. Default: all the field of view
addParameter(ip, 'StabilisationViewer', 1, @isnumeric)              %Stabilisation viewer. Default: yes

parse(ip, filename, varargin{:})
NameStabilisedVideo = ip.Results.NameStabilisedVideo;
algorithms = ip.Results.algorithms;
PercentualStrongestFeatures = ip.Results.PercentualStrongestFeatures;
TransformType = ip.Results.TransformType;
ROISelection = ip.Results.ROISelection;
StabilisationViewer = ip.Results.StabilisationViewer;

%% Reading Frames from a Video File
filename = num2str(filename);
VideoFrameReader = vision.VideoFileReader(filename, 'ImageColorSpace', 'Intensity');

obj = VideoReader(filename);
FPS = obj.FrameRate;
NFrames = floor(obj.Duration*FPS);

%% Processing all frames in the Video
reset(VideoFrameReader);
VideoPlayer = vision.VideoPlayer;   %Creating video viewer

% Initialising...
imgA = step(VideoFrameReader);      %First frame taken as the reference image
imgB = imgA;
imgBp = imgA;
correctedMean = imgA;

% Creating Video file to write stabilised frames
if length(algorithms) == 1
   NameStabilisedVideo = strcat(NameStabilisedVideo,'_',algorithms{1});
elseif length(algorithms) == 2
   NameStabilisedVideo = strcat(NameStabilisedVideo,'_',algorithms{1},'_',algorithms{2});
else
    error('ERROR!!! feature detection algorithms unclear')
end

v = VideoWriter(NameStabilisedVideo); %Name of the Stabilised Video
v.FrameRate = FPS;                    %Setting the original Frame Rate on the stabilised video
open(v)

writeVideo(v,imgA)                    %Saving reference image

% ROI Selection
if length(ROISelection) == 4
    ROI = ROISelection;
elseif length(ROISelection) == 1 && ROISelection == 1
    figure(1); imshow(imgA)            %Drawing a ROI
    r1 = drawrectangle('Label','ROI','Color',[1 0 0]);
    ROI = r1.Position;
    close;
elseif length(ROISelection) == 1 && ROISelection == 0
    ROI = NaN;
else
    error('ERROR!!! ROI Decision unclear')
end

% Analysis
Cont1 = 0; %Counting the frame number to filter geometric transformation
while ~isDone(VideoFrameReader)
    %% Reading a new frame
    imgAp = imgBp;
    imgB = step(VideoFrameReader);
    
    %% Estimating Geometric Transformations
    if length(ROISelection) == 4
        H = GeometricTransformStabilisation(imcrop(imgA,ROI),imcrop(imgB,ROI),algorithms,PercentualStrongestFeatures,TransformType);
    elseif length(ROISelection) == 1 && ROISelection == 1
        H = GeometricTransformStabilisation(imcrop(imgA,ROI),imcrop(imgB,ROI),algorithms,PercentualStrongestFeatures,TransformType);
    elseif length(ROISelection) == 1 && ROISelection == 0
        H = GeometricTransformStabilisation(imgA,imgB,algorithms,PercentualStrongestFeatures,TransformType);
    end
    
    %% Applying Geometric Transformations
    imgBp = imwarp(imgB,H,'OutputView',imref2d(size(imgA)));
    
    %% Saving Stabilised Frames
    writeVideo(v,imgBp)
    
    %% Displaying as color composite with last corrected frame
    if StabilisationViewer == 1
        step(VideoPlayer, imfuse(imgAp,imgBp,'ColorChannels','green-magenta'));
    elseif StabilisationViewer == 0
        Cont1 = Cont1 +1; clc
        ContProgress = 100*Cont1/NFrames
    else
        error('ERROR!!! Stabilisation Viewer Decision unclear')
    end
    correctedMean = correctedMean + imgBp;
end
close(v)

%% Realising Memory
release(VideoFrameReader);
release(VideoPlayer);
end


function H = GeometricTransformStabilisation(imgA,imgB,FeaturesDetection,StrongestFeatures_Decision,TransformType)
MultipleFeaturesDetection = size(FeaturesDetection);
SizeMFD = MultipleFeaturesDetection(1,2);

for NSizeMFD = 1:SizeMFD
    %% Features Detection
    DetAlgorithm = FeaturesDetection{1,NSizeMFD};
    FDA(1) = isequal(FeaturesDetection{1,NSizeMFD},char({'FAST'}));
    if FDA(1) == 1
        if  0 < StrongestFeatures_Decision && StrongestFeatures_Decision <= 100
            pointsA.(DetAlgorithm) = detectFASTFeatures(imgA);
            pointsB.(DetAlgorithm) = detectFASTFeatures(imgB);
            
            strongestPointsA.(DetAlgorithm) = selectStrongest(pointsA.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsA.(DetAlgorithm).Count));
            strongestPointsB.(DetAlgorithm) = selectStrongest(pointsB.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsB.(DetAlgorithm).Count));
            pointsA.(DetAlgorithm) = selectUniform(strongestPointsA.(DetAlgorithm),floor(0.5*strongestPointsA.(DetAlgorithm).Count),size(imgA));
            pointsB.(DetAlgorithm) = selectUniform(strongestPointsB.(DetAlgorithm),floor(0.5*strongestPointsB.(DetAlgorithm).Count),size(imgB));
        else
            error('ERROR!!! Percentage of strongest and uniform features out of limit (give a value between ]0, 100])')
        end
    end
    
    FDA(2) = isequal(FeaturesDetection{1,NSizeMFD},char({'MINEIGEN'}));
    if FDA(2) == 1
        if  0 < StrongestFeatures_Decision && StrongestFeatures_Decision <= 100
            pointsA.(DetAlgorithm) = detectMinEigenFeatures(imgA); % MinEigenFeatures: : Single-scale detection (point tracking with little or no scale change).
            pointsB.(DetAlgorithm) = detectMinEigenFeatures(imgB);
            
            strongestPointsA.(DetAlgorithm) = selectStrongest(pointsA.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsA.(DetAlgorithm).Count));
            strongestPointsB.(DetAlgorithm) = selectStrongest(pointsB.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsB.(DetAlgorithm).Count));
            pointsA.(DetAlgorithm) = selectUniform(strongestPointsA.(DetAlgorithm),floor(0.5*strongestPointsA.(DetAlgorithm).Count),size(imgA));
            pointsB.(DetAlgorithm) = selectUniform(strongestPointsB.(DetAlgorithm),floor(0.5*strongestPointsB.(DetAlgorithm).Count),size(imgB));
        else
            error('ERROR!!! Percentage of strongest and uniform features out of limit (give a value between ]0, 100])')
        end
    end
    
    FDA(3) = isequal(FeaturesDetection{1,NSizeMFD},char({'HARRIS'}));
    if FDA(3) == 1
        if  0 < StrongestFeatures_Decision && StrongestFeatures_Decision <= 100
            pointsA.(DetAlgorithm) = detectHarrisFeatures(imgA); % Harris: Single-scale detection (point tracking with little or no scale change).
            pointsB.(DetAlgorithm) = detectHarrisFeatures(imgB);
            
            strongestPointsA.(DetAlgorithm) = selectStrongest(pointsA.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsA.(DetAlgorithm).Count));
            strongestPointsB.(DetAlgorithm) = selectStrongest(pointsB.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsB.(DetAlgorithm).Count));
            pointsA.(DetAlgorithm) = selectUniform(strongestPointsA.(DetAlgorithm),floor(0.5*strongestPointsA.(DetAlgorithm).Count),size(imgA));
            pointsB.(DetAlgorithm) = selectUniform(strongestPointsB.(DetAlgorithm),floor(0.5*strongestPointsB.(DetAlgorithm).Count),size(imgB));
        else
            error('ERROR!!! Percentage of strongest and uniform features out of limit (give a value between ]0, 100])')
        end
    end
    
    FDA(4) = isequal(FeaturesDetection{1,NSizeMFD},char({'SURF'}));
    if FDA(4) == 1
        if  0 < StrongestFeatures_Decision && StrongestFeatures_Decision <= 100
            pointsA.(DetAlgorithm) = detectSURFFeatures(imgA); % SURF: Multiscale detection (Blobs - Object detection and image registration with scale and rotation changes)
            pointsB.(DetAlgorithm) = detectSURFFeatures(imgB);
            
            strongestPointsA.(DetAlgorithm) = selectStrongest(pointsA.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsA.(DetAlgorithm).Count));
            strongestPointsB.(DetAlgorithm) = selectStrongest(pointsB.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsB.(DetAlgorithm).Count));
            pointsA.(DetAlgorithm) = selectUniform(strongestPointsA.(DetAlgorithm),floor(0.5*strongestPointsA.(DetAlgorithm).Count),size(imgA));
            pointsB.(DetAlgorithm) = selectUniform(strongestPointsB.(DetAlgorithm),floor(0.5*strongestPointsB.(DetAlgorithm).Count),size(imgB));
        else
            error('ERROR!!! Percentage of strongest and uniform features out of limit (give a value between ]0, 100])')
        end
    end
    
    FDA(5) = isequal(FeaturesDetection{1,NSizeMFD},char({'KAZE'}));
    if FDA(5) == 1
        if  0 < StrongestFeatures_Decision && StrongestFeatures_Decision <= 100
            pointsA.(DetAlgorithm) = detectKAZEFeatures(imgA); % KAZE: Multi-scale blob features - Reduced blurring of object boundaries (point tracking with changes in scale in rotation)
            pointsB.(DetAlgorithm) = detectKAZEFeatures(imgB);
            
            strongestPointsA.(DetAlgorithm) = selectStrongest(pointsA.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsA.(DetAlgorithm).Count));
            strongestPointsB.(DetAlgorithm) = selectStrongest(pointsB.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsB.(DetAlgorithm).Count));
            pointsA.(DetAlgorithm) = selectUniform(strongestPointsA.(DetAlgorithm),floor(0.5*strongestPointsA.(DetAlgorithm).Count),size(imgA));
            pointsB.(DetAlgorithm) = selectUniform(strongestPointsB.(DetAlgorithm),floor(0.5*strongestPointsB.(DetAlgorithm).Count),size(imgB));
        else
            error('ERROR!!! Percentage of strongest and uniform features out of limit (give a value between ]0, 100])')
        end
    end
    
    FDA(6) = isequal(FeaturesDetection{1,NSizeMFD},char({'BRISK'}));
    if FDA(6) == 1
        if  0 < StrongestFeatures_Decision && StrongestFeatures_Decision <= 100
            pointsA.(DetAlgorithm) = detectBRISKFeatures(imgA); % BRISK: Multi-scale detection (point tracking with changes in scale in rotation)
            pointsB.(DetAlgorithm) = detectBRISKFeatures(imgB);
            
            strongestPointsA.(DetAlgorithm) = selectStrongest(pointsA.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsA.(DetAlgorithm).Count));
            strongestPointsB.(DetAlgorithm) = selectStrongest(pointsB.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsB.(DetAlgorithm).Count));
            pointsA.(DetAlgorithm) = selectUniform(strongestPointsA.(DetAlgorithm),floor(0.5*strongestPointsA.(DetAlgorithm).Count),size(imgA));
            pointsB.(DetAlgorithm) = selectUniform(strongestPointsB.(DetAlgorithm),floor(0.5*strongestPointsB.(DetAlgorithm).Count),size(imgB));
        else
            error('ERROR!!! Percentage of strongest and uniform features out of limit (give a value between ]0, 100])')
        end
    end
    
    FDA(7) = isequal(FeaturesDetection{1,NSizeMFD},char({'ORB'}));
    if FDA(7) == 1
        if  0 < StrongestFeatures_Decision && StrongestFeatures_Decision <= 100
            pointsA.(DetAlgorithm) = detectORBFeatures(imgA); % ORB: Multi-scale detection (point tracking with changes in scale in rotation)
            pointsB.(DetAlgorithm) = detectORBFeatures(imgB);
            
            strongestPointsA.(DetAlgorithm) = selectStrongest(pointsA.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsA.(DetAlgorithm).Count));
            strongestPointsB.(DetAlgorithm) = selectStrongest(pointsB.(DetAlgorithm),floor((StrongestFeatures_Decision/100)*pointsB.(DetAlgorithm).Count));
            pointsA.(DetAlgorithm) = selectUniform(strongestPointsA.(DetAlgorithm),floor(0.5*strongestPointsA.(DetAlgorithm).Count),size(imgA));
            pointsB.(DetAlgorithm) = selectUniform(strongestPointsB.(DetAlgorithm),floor(0.5*strongestPointsB.(DetAlgorithm).Count),size(imgB));
        else
            error('ERROR!!! Percentage of strongest and uniform features out of limit (give a value between ]0, 100])')
        end
    end
    
    if sum(FDA) == 0
        error('ERROR!!! No feature detection algorithm provided')
    end
    
    %% Feature correspondences
    % Extracting features
    [featuresA.(DetAlgorithm), ValidPointsA.(DetAlgorithm)] = extractFeatures(imgA, pointsA.(DetAlgorithm),'Method','Auto');
    [featuresB.(DetAlgorithm), ValidPointsB.(DetAlgorithm)] = extractFeatures(imgB, pointsB.(DetAlgorithm),'Method','Auto');
    
    % Matching features
    indexPairs.(DetAlgorithm) = matchFeatures(featuresA.(DetAlgorithm), featuresB.(DetAlgorithm));
    MatchedPointsA.(DetAlgorithm) = ValidPointsA.(DetAlgorithm)(indexPairs.(DetAlgorithm)(:, 1));
    MatchedPointsB.(DetAlgorithm) = ValidPointsB.(DetAlgorithm)(indexPairs.(DetAlgorithm)(:, 2));
end

if NSizeMFD == 1
    %% Use MSAC algorithm to compute the transformation
    H = estimateGeometricTransform(MatchedPointsB.(DetAlgorithm), MatchedPointsA.(DetAlgorithm), TransformType);
    
elseif NSizeMFD == 2
    matchedOriginalXY  = [MatchedPointsA.(char(FeaturesDetection{1,1})).Location; MatchedPointsA.(char(FeaturesDetection{1,2})).Location];
    matchedDistortedXY = [MatchedPointsB.(char(FeaturesDetection{1,1})).Location; MatchedPointsB.(char(FeaturesDetection{1,2})).Location];
    
    %% Use MSAC algorithm to compute the transformation
    H = estimateGeometricTransform(matchedDistortedXY, matchedOriginalXY, TransformType);
else
    error('ERROR!!! Maximum supported feature detection algorithms: 2')
end
end