function [valcube cbound nbound] = ContourTrack6(pre, post, tstr, tracksFinal, maxDiam, minD, minF, CFP, YFP, RFP, diagnostic)
% Sample
% valcube = ContourTrack6(pre, post, tstr, tracksFinal, 100, 1, .1, 933, 322, 125,0)
% Diagnostic = 1; Plot tracked nuclei with labels 
% Diagnostic = 2; Plot tracked nuclear and cytosolic bound


%Default Values for contour filtering
maxDiameter=maxDiam;
maxNucArea=round(pi*maxDiameter^2/4);
minDiameter=minD;
minNucArea=round(pi*minDiameter^2/4);
minFormfactor=minF;

% dilatesize=round((sqrt(sizeS*4)/pi)/2); %fix to scale with cell size
% dilatesize=round(sqrt(minNucArea/pi));
dilatesize=3;

%% Initialize video program if you want to take a video
if diagnostic==1
    vidObj = VideoWriter(['CDiagnostic1.avi']);           %Name
    open(vidObj);
    for f=1:length(tracksFinal)
        track=tracksFinal(f);
        tca=track.tracksCoordAmpCG;
        soe=track.seqOfEvents;
        xInd=1:8:length(tca);
        yInd=2:8:length(tca);
        startFrame=soe(1,1);
        timeRange=startFrame:(length(xInd)+startFrame-1);
        xCoord(f,timeRange)=tca(1,xInd);
        yCoord(f,timeRange)=tca(1,yInd);
    end
end

%% Get track coordinates at time from tracksFinal in a Matrix: xC and yC
for i=1:length(tracksFinal)
    %Get track coordinates in a matrix from tracksFinal.trackCoordAmpCG
    for j=1:length(tracksFinal(i).tracksCoordAmpCG)/8 - 1
        start = tracksFinal(i).seqOfEvents(1);
        %get x ((1+8n), row 1)
        xC(j+start,i) = tracksFinal(i).tracksCoordAmpCG(1+8*j);
        %get y ((2+8n), row 1)
        yC(j+start,i) = tracksFinal(i).tracksCoordAmpCG(2+8*j);
    end
end

xC(xC==0)=NaN;
yC(yC==0)=NaN;

xC2=xC(1:end,1:end);
yC2=yC(1:end,1:end);

% find NaN and fill with interpolated values
for i=1:size(xC,2)
    data=xC2(:,i);
    % If column of NaNs
    if isnan(max(data))==1
        xC(:,i)=NaN;
    else
        nanData = isnan(data);
        index = 1:numel(data);
        data2=data;
        % If only one value to interpolate
        if sum(~isnan(data2))==1
            xC(:,i)=data2;
            % Interpolate
        else
            data2(nanData) = interp1(index(~nanData), data(~nanData), index(nanData));
            xC(:,i)=data2;
        end
    end
end

for i=1:size(yC,2)
    data=yC2(:,i);
    % If column of NaNs
    if isnan(max(data))==1
        yC(:,i)=NaN;
    else
        nanData = isnan(data);
        index = 1:numel(data);
        data2=data;
        % If only one value to interpolate.
        if sum(~isnan(data2))==1
            yC(:,i)=data2;
            % Interpolate
        else
            data2(nanData) = interp1(index(~nanData), data(~nanData), index(nanData));
            yC(:,i)=data2;
        end
    end
end

%% Loop Structure to get intensities out of the tracks
for t=1:length(tstr)
% for t=289
    t
    % Read image
    filenameCFP=[pre{1} tstr{t} post{1}];
    filenameYFP=[pre{2} tstr{t} post{2}];
    filenameRFP=[pre{3} tstr{t} post{3}];
    im=imread(filenameYFP);
    % Get x and y coordinates for current track from current time.
    currentspots=find(xC(t,:));
    if ~isempty(currentspots)
    currentX=floor(xC(t,currentspots));
    currentY=floor(yC(t,currentspots));
    
    [maxY maxX] = size(im);
    
    % Turn coords into indices.
    for j=1:length(currentX)
        if currentX(j)>maxX
            currentX(j)=NaN;
        end
    end
    for i=1:length(currentY)
        if currentY(i)>maxY
            currentY(i)=NaN;
        end
    end
    
    currentInd=sub2ind(size(im),currentY,currentX);
    %Get filtered index (nonzero value tracks)
    filInd=currentInd(currentInd>0);
    
    % Create contour matrix
    er=reshape(im,size(im,1)*size(im,2),1);                              	
    q1=double(quantile(er,0.05));                                           
    q9=double(quantile(er,0.98));
    thresholds=fliplr(linspace(q1,q9,30));
    ctm=zeros(size(im));
    k=1;
    
    % While loop to filter best contour slice for each coordinate. Start at
    % biggest contour and work way down. Break while loop when first contour 
    % that satisfies conditions is met.
    
    while ~isempty(filInd) && k<=length(thresholds)
%             tim=bwlabel(imfill(~(im>thresholds(k)),'holes'));
        tim=bwlabel(~(im>thresholds(k)));
        ct = imopen(tim,strel('disk',2));
        currentLabels=ct(filInd);                                           % find the labels of the regions in which current Indices fall
        ct=ismember(ct,currentLabels(currentLabels>1));                     % remove regions in which no indies fall; ct is now binary
        ct=bwlabel(ct);                                                     % convert ct back to a label matrix
        currentLabels=ct(filInd);                                           % get the labels for the current Indices
        S = regionprops(ct,'EquivDiameter','Area','Perimeter');
        nucArea=cat(1,S.Area);
        nucPerim=cat(1,S.Perimeter);
        nucEquiDiameter=cat(1,S.EquivDiameter);
        nucFormfactor=4*pi*nucArea./(nucPerim.^2);
        sizeScore = nucArea>minNucArea*.3 & nucArea<maxNucArea;
        shapeScore=nucFormfactor>minFormfactor;
        totalScore=sizeScore&shapeScore;
        scorenz=find(totalScore);
        % Keep contour for that satisfies condition for current coord.
        ctm(ismember(ct,scorenz))=1;
        % Remove found contours from list of indices to search for
        foundLabels=unique(ct(ismember(ct,scorenz)));
        filInd(ismember(currentLabels,foundLabels))=[];
        k=k+1;
    end
    
    currentInd=sub2ind(size(im),currentY,currentX);
    goodInd=find(currentInd>0);
    ctml=bwlabel(ctm);
    nlabels=(ctml(currentInd(goodInd)));
    
    CFPim = imread([pre{1} tstr{t} post{1}]);
    YFPim = imread([pre{2} tstr{t} post{2}]);
    if RFP~=0
    RFPim = imread([pre{3} tstr{t} post{3}]);
    end
    % if you want to plot the tracked cells with labels.
    if diagnostic==1
        h=figure(1); clf; hold on;
        imshow(im,[]); hold on;
                for p=1:length(tracksFinal)
                    plot(xCoord(p,t),yCoord(p,t),'c.');
                    text(xCoord(p,t),yCoord(p,t),num2str(p),'Color','Green','FontSize',10);
                end
        set(0,'CurrentFigure',1);
        frame = getframe;
        writeVideo(vidObj,frame);
        
    end
    % If you want to see the tracked cells nuclear and cytosolic bounds
    if diagnostic ==2
    figure(1); clf; hold on;
    imshow(im,[]); hold on;
    end
    for n=1:length(goodInd)
        % Get intensity values for each region.
        if nlabels(n)~=0
            nmask=ctml==nlabels(n);
            nb=bwboundaries(nmask);
            nbound(t,goodInd(n))=nb(1);
            regOnly=ctml==nlabels(n);
            se=strel('disk',dilatesize);
            cmask=imdilate(regOnly,se);
            cmask(nmask)=0;
            cb=bwboundaries(cmask);
            cbound(t,goodInd(n))=cb(1);
            % Get intensities for nuclear region
            valcube(goodInd(n),t,1)=mean(CFPim(nmask));
            valcube(goodInd(n),t,2)=mean(YFPim(nmask));
            if RFP~=0
            valcube(goodInd(n),t,3)=mean(RFPim(nmask));
            end
            %Get intensities for donut region
            valcube(goodInd(n),t,4)=mean(CFPim(cmask));
            valcube(goodInd(n),t,5)=mean(YFPim(cmask));
            if RFP~=0
            valcube(goodInd(n),t,6)=mean(RFPim(cmask));
            end
            if diagnostic==1 | 2
                %Plot dilation
                mycellc=cb(size(cb,1),:);
                if isempty(mycellc{1})
                    continue
                else
                    dim=mycellc{1};
                    dim2 = cat(2, dim(:,1),dim(:,2));
                    hold on, plot(dim2(:,2),dim2(:,1),'r-');
                end
                %Plot nuclear bound
                mycelln=nb(size(nb,1),:);
                if isempty(mycelln{1})
                    continue
                else
                    dim=mycelln{1};
                    dim2 = cat(2, dim(:,1),dim(:,2));
                    hold on, plot(dim2(:,2),dim2(:,1),'g-');
                end
            end
        else
            %valcube(goodInd(n,t,1)=valcube(goodInd(n-1),(t-1),1);
            valcube(goodInd(n),t,1)=NaN;
            valcube(goodInd(n),t,2)=NaN;
            valcube(goodInd(n),t,3)=NaN;
            %Get intensities for donut region
            valcube(goodInd(n),t,4)=NaN;
            valcube(goodInd(n),t,5)=NaN;
            valcube(goodInd(n),t,6)=NaN;
        end
    end
    else
    end
    
    
end
% Close video
if diagnostic==1
    close(vidObj);
end

%% Get the ratio intensities
% Get CFP/YFP intensity ratio for nuclear region
valcube(:,:,7)=(valcube(:,:,1))./(valcube(:,:,2));
% Get CFP/YFP intensity for donut region
valcube(:,:,8)=(valcube(:,:,4))./(valcube(:,:,5));
cbg=valcube(:,:,4)-CFP;
ybg=valcube(:,:,5)-YFP;
ratio=cbg./ybg;
filter=cbg>0;
ratio2=ratio.*filter;
valcube(:,:,9)=ratio2;
% Get CFP/RFP intensity ratio for nuclear region
valcube(:,:,10)=(valcube(:,:,1))./(valcube(:,:,3));
% Get CFP/RFP intensity for donut region
if RFP~=0
valcube(:,:,11)=(valcube(:,:,4))./(valcube(:,:,6));
cbg=valcube(:,:,4)-CFP;
ybg=valcube(:,:,6)-RFP;
ratio3=cbg./ybg;
filter=cbg>0;
ratio4=ratio3.*filter;
valcube(:,:,12)=ratio4;
filenameCFP=[pre{1} tstr{1} post{1}];
end
save(['data_xy2' filenameCFP '.mat']);

