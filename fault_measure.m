% fault measuring

% This script calls necessary companion functions to read and automatically
% measure fault throws from grdtrack (GMT) outputs

%

% ***** RECORD *****
% 24 Apr 2018 | first cut | T.Morrow
%
%
% ***** ****** *****

% Input format should be series of cross-profiles from grdtrack,
% instructions below:

%
%   After generating a .xy file for one fault, I use the following command
%       grdtrack FT1.xy -GFAULT.grd -C4000e/200e/200e > FTXprof1.txt
%                                     ^1    ^2   ^3
%           arguments 1, 2, 3 are cross profile length, along-cross
%           spacing, and along-fault interval respectively
%   Then, concatenate all the files together
%       cat FTXprof* >> allfts.txt
%
%   Import this result into Matlab as a single matrix (will have NaNs
%   between each segment)
%
%   DEPENDENCIES
%   GMT used to make input
%   grdread2 used in demo
%
%
%   OUTPUTS
%   Heave and Throw arrays, with one measurement (m) per cross-section.
%   NOTE: this method yields good measurements for throw but not so great
%   for heave. Run demo case (see below) to see sample measurement (pink
%   line)
%
%   DEMO for one cross-section available at bottom of code

%% CODE

% naming conventions
%   allfts    imported fault numeric matrix

close all

% index cross profiles
I=find(isnan(allfts(:,1)));

% parse faults to cell array
for qq=1:length(I)
    if qq==1
        FTXs{qq}=allfts(1:I(qq)-1,:);
    elseif qq==length(I)
        FTXs{qq}=allfts(I(qq-1)+1:end,:);
    else
        FTXs{qq}=allfts(I(qq-1)+1:I(qq)-1,:);
    end
end

% measure each fault cross section
for rr=1:length(I)
    XP=FTXs{rr};

    % matrix of XYZ pts defining one bathymetric profile across a fault
    XFT = [XP(:,3) zeros(size(XP(:,1))) XP(:,5)];

    XFT(:,2:4)=XFT;
    XFT(:,1) = 1:length(XFT(:,2));

    XFT = double(XFT);

    % bathymetry resolution
    rez = 60;
    %plot3(XFT(:,2),XFT(:,3),XFT(:,4),'b');
    %hold all

    % minimum and maximum dip of faults
    %DIPMN = minDIP;
    %DIPMX = 90;

    % standard deviation window size
    w=40;

    %% Calculate gradient along cross profile
    GRDFT = horzcat(XFT(:,1:3), abs(gradient(XFT(:,4))));
    %plot3(GRDFT(:,2),GRDFT(:,3),GRDFT(:,4)+mean(XFT(:,4)),'r') % -plot test for gradient- %

    %% Find largest peak in standard deviation [not currently using]
    T = XFT(:,4);
    stdT=T.*0;

    for i=1+w/2:length(stdT)-w/2-1
        im=i-w/2:i+w/2;
        stdT(i)=nanstd(T(im));
    end

    % find peaks
    STD = horzcat(XFT(:,1:3), stdT);
    %plot3(STD(:,2),STD(:,3),STD(:,4)+mean(XFT(:,4)),'k');
    STD_PKS = findpeaks(STD(:,4));

    % find peak locations
    for ii = 1:length(STD_PKS)
         if (ii == 1);
            STD_LOC_old = STD(find(STD(:,4)==STD_PKS(ii,1)),1:4);
        else
            STD_LOC_new = STD(find(STD(:,4)==STD_PKS(ii,1)),1:4);
            STD_LOC_old = vertcat(STD_LOC_old, STD_LOC_new);
         end
    end

    STD_LOC = flipud(sortrows(STD_LOC_old,4));
    %plot3(STD_LOC(:,2),STD_LOC(:,3),ones(length(STD_LOC),1)*mean(XFT(:,4)),'^k'); % -plot test for std peaks- %
    MXSTD = STD_LOC(1,1:3);

    %% Find peak-based bounds
    GRD_BND = horzcat(GRDFT(:,1:3), -GRDFT(:,4));
    %plot3(GRD_BND(:,2),GRD_BND(:,3),GRD_BND(:,4)+mean(XFT(:,4)),'r') % -plot test for inv gradient- %

    % find peaks of inverted gradient
    BND_PKS = findpeaks(GRD_BND(:,4));
    for ii = 1:length(BND_PKS);
        if (ii == 1)
            BNDPK_LOC_old = GRD_BND(find(GRD_BND(:,4)==BND_PKS(ii,1)),1:3);
        else
            BNDPK_LOC_new = GRD_BND(find(GRD_BND(:,4)==BND_PKS(ii,1)),1:3);
            BNDPK_LOC_old = vertcat(BNDPK_LOC_old, BNDPK_LOC_new);
        end
    end
    BNDPK_LOC = BNDPK_LOC_old;

    BNDPK_LOC = sortrows(BNDPK_LOC,1);
    %plot3(BNDPK_LOC(:,2),BNDPK_LOC(:,3),ones(length(BNDPK_LOC(:,1)),1)*mean(XFT(:,4)),'*b'); % -plot test for inv gradient peaks- %

    % flag and concatenate peaks on gradient and inv gradient
    MXSTD(:,4) = ones(length(MXSTD(:,1)),1);
    BNDPK_LOC(:,4) = zeros(length(BNDPK_LOC(:,1)),1);

    ALLPKS = vertcat(MXSTD, BNDPK_LOC);
    ALLPKS = sortrows(ALLPKS,1);

    % locate bounding points for each fault (nearest peak in inv gradient)
    FTind = find(ALLPKS(:,4));

    
    %% Take each type of peak and find greatest vertical offset combination
    
    for ii=1:length(BNDPK_LOC(:,1))
        RAT_LOC(ii,:)=XFT(BNDPK_LOC(ii,1),:);
    end

    RAT=[];
    pp=1;
    for ii=1:length(RAT_LOC(:,1))
        for jj=1:length(RAT_LOC(:,1))
            if ii==jj
                continue
            else
                %RAT(pp,:)=[ii jj (RAT_LOC(ii,2)-RAT_LOC(jj,2))/(RAT_LOC(ii,4)-RAT_LOC(jj,4))^2];
                RAT(pp,:)=[ii jj 1/(RAT_LOC(ii,4)-RAT_LOC(jj,4))];
                pp=pp+1;
            end
        end
    end

    RAT=sortrows(abs(RAT),3);

    XFT_pts=[RAT_LOC(RAT(1,1),:); RAT_LOC(RAT(1,2),:)];
    
    %plot3(XFT_pts(:,2),[0;0],XFT_pts(:,4),'m-o')
    
    %% Measure faults
    
    HV = abs(XFT_pts(1,2)-XFT_pts(2,2));
    THRO = abs(XFT_pts(1,4)-XFT_pts(2,4));
    %for ll=1:length(BND);
    %    HV{ll} = sqrt((BND{ll}(1,2)-BND{ll}(2,2))^2 + (BND{ll}(1,3)-BND{ll}(2,3))^2);
    %    THRO{ll} = abs(XFT(BND{ll}(1,1),4)-XFT(BND{ll}(2,1),4));
    %end

    %HV_TOT = sum(vertcat(HV{:}));
    %THRO_TOT = sum(vertcat(THRO{:}));
    %display(['Total heave ',num2str(HV_TOT),'m']);
    %display(['Total throw ',num2str(THRO_TOT),'m']);
    
    Heave{rr}=HV;
    Throw{rr}=THRO;
    
    display(['Heave ',num2str(HV),'m']);
    display(['Throw ',num2str(THRO),'m']);
    
    clear BND STD STD_LOC STD_LOC_new STD_LOC_old THRO HV XFT MXSTD BNDPK_LOC BNDPK_LOC_new BNDPK_LOC_old RAT RAT_LOC STD_PKS XFT_pts GRD_BND GRDFT ALLPKS
end

%% DEMO BELOW REMOVE LAST LINE OF CODE AND LINE BELOW TO ACTIVATE
%{

%% CODE

% naming conventions
%   allfts    imported fault numeric matrix
%   X,Y,Z     bathymetry read in from grdread2

% load sample data set
load sample_imported_data.mat

% read in bathymetry (replace with your grid file)
[X Y Z]=grdread2('FAULT.grd');

close all

    % plot
    subplot(121)
    contourf(X,Y,Z,20);
    axis equal
    colorbar
    hold on
    
%figure

% index cross profiles
I=find(isnan(allfts(:,1)));

% parse faults to cell array
for qq=1:length(I)
    if qq==1
        FTXs{qq}=allfts(1:I(qq)-1,:);
    elseif qq==length(I)
        FTXs{qq}=allfts(I(qq-1)+1:end,:);
    else
        FTXs{qq}=allfts(I(qq-1)+1:I(qq)-1,:);
    end
end

for zz=1:length(FTXs)
    plot(FTXs{zz}(FTXs{zz}(:,3)==0,1),FTXs{zz}(FTXs{zz}(:,3)==0,2),'rx');
end

% measure one random fault cross section
for rr=round(rand(1)*length(I)) %1:length(I)
    XP=FTXs{rr};

    plot([XP(1,1) XP(end,1)],[XP(1,2) XP(end,2)],'g-*','LineWidth',3)
    
    subplot(122)
    hold on
    
    
    % matrix of XYZ pts defining one bathymetric profile across a fault
    XFT = [XP(:,3) zeros(size(XP(:,1))) XP(:,5)];

    XFT(:,2:4)=XFT;
    XFT(:,1) = 1:length(XFT(:,2));

    XFT = double(XFT);

    % bathymetry resolution
    rez = 60;
    plot(XFT(:,2),XFT(:,4),'b');
    %hold all

    % minimum and maximum dip of faults
    %DIPMN = minDIP;
    %DIPMX = 90;

    % standard deviation window size
    w=40;

    %% Calculate gradient along cross profile
    GRDFT = horzcat(XFT(:,1:3), abs(gradient(XFT(:,4))));
    %plot3(GRDFT(:,2),GRDFT(:,3),GRDFT(:,4)+mean(XFT(:,4)),'r') % -plot test for gradient- %

    %% Find largest peak in standard deviation [not currently using]
    T = XFT(:,4);
    stdT=T.*0;

    for i=1+w/2:length(stdT)-w/2-1
        im=i-w/2:i+w/2;
        stdT(i)=nanstd(T(im));
    end

    % find peaks
    STD = horzcat(XFT(:,1:3), stdT);
    %plot3(STD(:,2),STD(:,3),STD(:,4)+mean(XFT(:,4)),'k');
    STD_PKS = findpeaks(STD(:,4));

    % find peak locations
    for ii = 1:length(STD_PKS)
         if (ii == 1);
            STD_LOC_old = STD(find(STD(:,4)==STD_PKS(ii,1)),1:4);
        else
            STD_LOC_new = STD(find(STD(:,4)==STD_PKS(ii,1)),1:4);
            STD_LOC_old = vertcat(STD_LOC_old, STD_LOC_new);
         end
    end

    STD_LOC = flipud(sortrows(STD_LOC_old,4));
    %plot3(STD_LOC(:,2),STD_LOC(:,3),ones(length(STD_LOC),1)*mean(XFT(:,4)),'^k'); % -plot test for std peaks- %
    MXSTD = STD_LOC(1,1:3);

    %% Find peak-based bounds
    GRD_BND = horzcat(GRDFT(:,1:3), -GRDFT(:,4));
    %plot3(GRD_BND(:,2),GRD_BND(:,3),GRD_BND(:,4)+mean(XFT(:,4)),'r') % -plot test for inv gradient- %

    % find peaks of inverted gradient
    BND_PKS = findpeaks(GRD_BND(:,4));
    for ii = 1:length(BND_PKS);
        if (ii == 1)
            BNDPK_LOC_old = GRD_BND(find(GRD_BND(:,4)==BND_PKS(ii,1)),1:3);
        else
            BNDPK_LOC_new = GRD_BND(find(GRD_BND(:,4)==BND_PKS(ii,1)),1:3);
            BNDPK_LOC_old = vertcat(BNDPK_LOC_old, BNDPK_LOC_new);
        end
    end
    BNDPK_LOC = BNDPK_LOC_old;

    BNDPK_LOC = sortrows(BNDPK_LOC,1);
    %plot3(BNDPK_LOC(:,2),BNDPK_LOC(:,3),ones(length(BNDPK_LOC(:,1)),1)*mean(XFT(:,4)),'*b'); % -plot test for inv gradient peaks- %

    % flag and concatenate peaks on gradient and inv gradient
    MXSTD(:,4) = ones(length(MXSTD(:,1)),1);
    BNDPK_LOC(:,4) = zeros(length(BNDPK_LOC(:,1)),1);

    ALLPKS = vertcat(MXSTD, BNDPK_LOC);
    ALLPKS = sortrows(ALLPKS,1);

    % locate bounding points for each fault (nearest peak in inv gradient)
    FTind = find(ALLPKS(:,4));

    
    %% Take each type of peak and find greatest vertical offset combination
    
    for ii=1:length(BNDPK_LOC(:,1))
        RAT_LOC(ii,:)=XFT(BNDPK_LOC(ii,1),:);
    end

    RAT=[];
    pp=1;
    for ii=1:length(RAT_LOC(:,1))
        for jj=1:length(RAT_LOC(:,1))
            if ii==jj
                continue
            else
                %RAT(pp,:)=[ii jj (RAT_LOC(ii,2)-RAT_LOC(jj,2))/(RAT_LOC(ii,4)-RAT_LOC(jj,4))^2];
                RAT(pp,:)=[ii jj 1/(RAT_LOC(ii,4)-RAT_LOC(jj,4))];
                pp=pp+1;
            end
        end
    end

    RAT=sortrows(abs(RAT),3);

    XFT_pts=[RAT_LOC(RAT(1,1),:); RAT_LOC(RAT(1,2),:)];
    
    plot(XFT_pts(:,2),XFT_pts(:,4),'m-o')
    
    %% Measure faults
    
    HV = abs(XFT_pts(1,2)-XFT_pts(2,2));
    THRO = abs(XFT_pts(1,4)-XFT_pts(2,4));
    %for ll=1:length(BND);
    %    HV{ll} = sqrt((BND{ll}(1,2)-BND{ll}(2,2))^2 + (BND{ll}(1,3)-BND{ll}(2,3))^2);
    %    THRO{ll} = abs(XFT(BND{ll}(1,1),4)-XFT(BND{ll}(2,1),4));
    %end

    %HV_TOT = sum(vertcat(HV{:}));
    %THRO_TOT = sum(vertcat(THRO{:}));
    %display(['Total heave ',num2str(HV_TOT),'m']);
    %display(['Total throw ',num2str(THRO_TOT),'m']);
    
    Heave{rr}=HV;
    Throw{rr}=THRO;
    
    display(['Heave ',num2str(HV),'m']);
    display(['Throw ',num2str(THRO),'m']);
    
    clear BND STD STD_LOC STD_LOC_new STD_LOC_old THRO HV XFT MXSTD BNDPK_LOC BNDPK_LOC_new BNDPK_LOC_old RAT RAT_LOC STD_PKS XFT_pts GRD_BND GRDFT ALLPKS
end

%}
