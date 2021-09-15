% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This script enables the mainP2.m file to be a parametric script to run multiple simulations

function [Energy, PRF, primaryLens, eMF, eENF, rIF_BW, rMF, anodeDarkCurrent, anodeRadiantSens, vsm, LidarRatio, boundaryLayer, opertaionTime] = parameters(code)

if code == 0
    Energy = 160;
    PRF = 20;
    primaryLens = 0.2031;
    eMF = 150;
    eENF = 4.5;
    rIF_BW = 3;
    rMF = 3*10^6;
    anodeDarkCurrent = 1;
    anodeRadiantSens = 3*10^4;
    vsm = 39.12;
    LidarRatio = 25;
    boundaryLayer = 3;
    opertaionTime = 1; %nightTime    
end

if code == 1
    Energy = 300;
    PRF = 40;
    primaryLens = 0.5;
    eMF = 150;
    eENF = 4.5;
    rIF_BW = 0.5;
    rMF = 3*10^6;
    anodeDarkCurrent = 1;
    anodeRadiantSens = 3*10^4;
    vsm = 3.912;
    LidarRatio = 25;
    boundaryLayer = 3;
    opertaionTime = 2; %dayTime    
end

if code == 2
    Energy = 100;
    PRF = 100;
    primaryLens = 0.4;
    eMF = 400;
    eENF = 4.5;
    rIF_BW = 3; %3
    rMF = 3*10^6;
    anodeDarkCurrent = 1;
    anodeRadiantSens = 6*10^4;
    vsm = 39.12;
    LidarRatio = 25;
    boundaryLayer = 3;
    opertaionTime = 2; %dayTime    
end

if code == 3
    Energy = 150;
    PRF = 30;
    primaryLens = 0.3;
    eMF = 150;
    eENF = 4.5;
    rIF_BW = 3;
    rMF = 1*10^6;
    anodeDarkCurrent = 1;
    anodeRadiantSens = 3*10^4;
    vsm = 3.912;
    LidarRatio = 25;
    boundaryLayer = 3;
    opertaionTime = 1; %nightTime    
end

if code == 5
    Energy = 150;
    PRF = 30;
    primaryLens = 0.3;
    eMF = 150;
    eENF = 4.5;
    rIF_BW = 3;
    rMF = 1*10^6;
    anodeDarkCurrent = 1;
    anodeRadiantSens = 3*10^4;
    vsm = 3.912;
    LidarRatio = 25;
    boundaryLayer = 3;
    opertaionTime = 1; %nighTime
end

if code == 6
    Energy = 100;
    PRF = 100;
    primaryLens = 0.4;
    eMF = 150;
    eENF = 4.5;
    rIF_BW = 0.5;
    rMF = 3*10^6;
    anodeDarkCurrent = 1;
    anodeRadiantSens = 3*10^4;
    vsm = 39.12;
    LidarRatio = 25;
    boundaryLayer = 2;
    opertaionTime = 1; %nighTime
end

end

