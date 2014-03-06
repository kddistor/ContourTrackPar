function [movieInfo,ecube,etlm]=qte4(namePre,namePost,tstr,trackChannel,maxDiam, minD, minF, diagnostic)
% Sample input
% movieInfo = qte4(pre,post,tstr,2,100,14,.6,1);
% Diagnostic = 1; Plot the nuclear bounds of identified cells

%% Get parameters for Nuclei Detection
maxDiameter=maxDiam;
maxNucArea=round(pi*maxDiameter^2/4);
minDiameter=minD;
minNucArea=round(pi*minDiameter^2/4);
minFormfactor=minF;
timepoints=length(tstr);
hy = fspecial('sobel');                                                    %Create a predefined 2d filter. "Sobel horizontal edge-emphasizing filter" Can use Canny
hx = hy';
%% Initialize diagnostic
if diagnostic==1
filename=[namePre{trackChannel} tstr{1} namePost{trackChannel}];
vidObj = VideoWriter([filename 'Diagnostic1.avi']);           %Name
open(vidObj);
end

%% Loop Structure to detect cell nuclei positions
for t=1:length(tstr)
    t
    filename=[namePre{trackChannel} tstr{t} namePost{trackChannel}];        %Open Image.
    im=imread(filename);                                                    %Reads Image.
%     background = imopen(im,strel('disk',7));                               %Create background structuring element
%     im = imsubtract(im,background);                                                  %Subtract background
%     gaussianFilter = fspecial('gaussian', [3, 3], 5);                      %Gaussian filter to filter noise              
    
%Added 2/25 3:10pm MPargett ----------------------------------------
    %Gaussian filter (smoothes noise)
    gaussianFilter = fspecial('gaussian', [3, 3], 10);
    im2 = imfilter(im, gaussianFilter, 'replicate');
    
    % Remove background
    nbins = numel(unique(im2(:)))*sqrt(2); %Get number of unique values
    %Robustly identify first peak of histogram
    for s = 1:100
        %Iterate nbins
        nbins = ceil(nbins./sqrt(2));
        %Get histogram of image
        [hh, xx] = hist( double(im2(:)), nbins );
        %Check for empty bins in front tail
        fi = find(hh,1,'first');
        if ~any( hh(fi:fi+floor(nbins/10)) == 0 ); break; end
        %   Note: heurisitc value used for 'low' region
    end
    %Smooth histogram to exclude small variations
    hh = smooth(hh, 5, 'moving');
    hh = [0; hh]; %#ok<AGROW> %Pad in case first bin is a peak
    %Find extrema (where derivative sign changes)
    extma = find( sign(hh(3:end) - hh(2:end-1)) ~=...
        sign(hh(2:end-1) - hh(1:end-2)) ) + 1;
    %Peaks only (exclude minima, by second derivative)
    pks = xx(extma(hh(extma+1) + hh(extma-1) - 2*hh(extma) < 0) - 1);
    %Remove background by subtracting first peak
    im2 = im2 - pks(1);  %Subtract first peak
%---------------------------------------------------------------------

    Iy = imfilter(double(im2), hy, 'replicate');                             %N-D filtering on 2D filtered image.
    Ix = imfilter(double(im2), hx, 'replicate');                             %Across x axis.
    e = sqrt(Ix.^2 + Iy.^2);                                                %Create vector of distances between corresponding points.
    er=reshape(e,size(e,1)*size(e,2),1);                                    %Reshape array according to distance by sizes of distance of both points.
    q1=quantile(er,0.05);                                                 	%Find the value at which 10% of points falls beneath.
    q9=quantile(er,0.9);                                                    %Find the value at which 90% of points falls beneath.
    thresholds=fliplr(linspace(q1,q9,10));                                  %Create a vector of from quantile values
    etlm=zeros(size(im2));                                                   %Initialize label matrix for each pixel
    for j=1:length(thresholds)
        et=e>thresholds(j);                                                 %Create an vector from distances that are greater than threshholds
%         eti=~et;                                                            %Invert image
%         etl=bwlabel(double(eti)./double(im));                               %Divide image by original to get brighter nuclei
        et = imdilate(et, strel('disk',1));                                %Dilate image
        et = imerode(et, strel('disk',2));                                 %Erode image to get rid of background noise.
%         et = imerode(et, strel('disk',2));                                 %Erode image to get rid of background noise.
        et = imdilate(et, strel('disk',1));                                %Dilate image
%         et = bwareaopen(et, minNucArea);                                  %Get rid of small obects
        etl=bwlabel(~et);                                                   %Invert image black/white
        S = regionprops(etl,'EquivDiameter','Area','Perimeter');            %Finds region properties of CC and stores into S.
        nucArea=cat(1,S.Area);                                              %Gets the Area of S.
        nucPerim=cat(1,S.Perimeter);                                        %Gets Perimeter of S.
        nucEquiDiameter=cat(1,S.EquivDiameter);                             %Gets the diameter of a circle with the same area as the region. Computed as sqrt(4*Area/pi).
        nucFormfactor=4*pi*nucArea./(nucPerim.^2);                          %Gets the shape value for threshholding.
        sizeScore = nucArea>minNucArea*0.3 & nucArea<maxNucArea;            %Returns 1 if both statements are true. Else returns 0.
        shapeScore=nucFormfactor>minFormfactor;                             %Returns 1 if statement is true. Else returns 0.
        totalScore=sizeScore&shapeScore;                                    %Returns 1 if both statements are true. Else returns 0.
        scorenz=find(totalScore);                                           %Score index
        for k=1:length(scorenz)                                             %Filter by score
            etlm(etl==scorenz(k))=1;                                        %Find indexes where score is good.
        end
    end
    C=regionprops(etlm>0,'Centroid','Area');                                %Get Region Stats.
    %Restructure all points to be piped to uTrack2.0 (scriptTrackGeneral)
    P1 = cat(1, C.Centroid);                               
    P2 = cat(1, C.Area);
    if ~isempty(P1)
        z=zeros(length(P1),1);
        P1(:,1);
        if size(P1,1)==1
            z=zeros(1,1);
        end
        xC=cat(2,P1(:,1),z);
        yC=cat(2,P1(:,2),z);
        am=cat(2,P2(:,1),z);
        s1.xCoord=xC;
        s1.yCoord=yC;
        s1.amp=am;
        
        mI(t)=s1;
        ecube(:,:,t)=etlm;
    end
%% Create movie if diagnostic is on
    if diagnostic==1
        h=figure(1); clf; hold on;
        imshow(im,[]); hold on;
        nmask=etlm;
        nb=bwboundaries(nmask);
        for c=1:size(nb,1)
            mycell=nb(c,:);
            if isempty(mycell{1})
                continue
            else
                dim=mycell{1};
                dim2 = cat(2, dim(:,1),dim(:,2));
                hold on, plot(dim2(:,2),dim2(:,1),'y-');
            end
        end
    set(0,'CurrentFigure',1);
    frame = getframe;
    writeVideo(vidObj,frame);
    end
end
if diagnostic==1
close(vidObj);
end
movieInfo=mI;