%%  Name: HLTVSG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT NOTE:
% This demo is applicable to the noise case given in the HLTVSG paper. 
% We have provided an experimental case of the ROSIS Pavia city center data case 4 and case6 in the HLTVSG_Demo.
% The default parameters are fit for the noise cases in the paper.
% It is worth noting that when the noise intensity becomes very high, there should be some subtle modification of the parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath(genpath('HLTVSG'));
addpath(genpath('lib'));
addpath(genpath('quality_assess'));
addpath(genpath('tensor_toolbox'));
addpath(genpath('TV_operator'));

load NPavia_case6.mat
[N_psnr, N_ssim, N_msam] = MSIQA3(OriData3*255, oriData3_noise*255);
r = 2; 
beta = 3.8; 
lambda = 2.9; 
tic;
HLTVSG_output = HLTVSG(oriData3_noise, beta,lambda, r);
HLTVSG_time = toc;
[HLTVSG_psnr, HLTVSG_ssim, HLTVSG_msam] = MSIQA3(OriData3*255, HLTVSG_output*255);