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


%% Initialising
clear all
close all
clc

%% Defining Variables
%Input:
filename = 'Belgrade_15frames.avi';     %Name of the video to be stabilised

%Configurable Settings:
NameStabilisedVideo = 'TestVideo';      %Name of the stabilised video. Default: 'StabilisedVideo + Feature detection algorithms'. 
algorithms = {'FAST'};                  %Features' detection algorithms (Max 2). Default: {'SURF'}.
PercentualStrongestFeatures = 35;       %Percentage value of the strongest detected features (Value between ]0, 100]). Default: 100.
TransformType = 'projective';           %Transform type. Default: 'similarity'.  
ROISelection = [1 1 1000 1000];         %ROI decision. Default: all the field of view.
StabilisationViewer = 0;                %Stabilisation viewer. Default: 1 (viewer)

%% Calling VISION (run it one by one to see the differences | comment and uncomment when necessary)
%Running VISION with default settings:
[NFrames,FPS,ROI] = VISION(filename)

% %Giving a different name to the stabilised video:
% [NFrames,FPS,ROI] = VISION(filename,'NameStabilisedVideo',NameStabilisedVideo)

% %Using a different feature detection algorithm:
% [NFrames,FPS,ROI] = VISION(filename,'algorithms',algorithms)

% %Using a different percentage of strongest features:
% [NFrames,FPS,ROI] = VISION(filename,'PercentualStrongestFeatures',PercentualStrongestFeatures)

% %Using a different transform type:
% [NFrames,FPS,ROI] = VISION(filename,'TransformType',TransformType)

% %Using a vector ROI:
% [NFrames,FPS,ROI] = VISION(filename,'ROISelection',ROISelection)

% %Using a drawned ROI:
% [NFrames,FPS,ROI] = VISION(filename,'ROISelection',1)

% %Running VISION without stabilisation viewer:
% [NFrames,FPS,ROI] = VISION(filename,'StabilisationViewer',StabilisationViewer)

