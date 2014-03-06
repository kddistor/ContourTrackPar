for h=7:9
    xystring= ['xy0', num2str(h)]
    [pre,post,tstr]=getNames(3,1,1,484,1,5,xystring,1,2,'2014-02-15-ekarfra','c1','c2','c3','c4','tif',[10 50]);
    parfor i=1:length(tstr)
        movieInfo{i} = qtetest(pre,post,tstr(i),2,100,14,.6,0);
    end
    movieInfo = cat(2, movieInfo{:});
    scriptTrackGeneral
    %Set the number of workers you want to use
    %You should indicate more workers than you have data elements (I don't think the code is robust to that)
    nWorkers = 9;
    
    %For your application, with tstr (this could be any data that you wanted to splice)
    nTime = numel(tstr);
    
    %Pre-define the cell containing the spliced data
    tstru = cell(nWorkers,1);
    %Calculate the max number of points to put in each cell
    nPer = ceil(nTime/nWorkers);
    
    for s = 1:nWorkers
        %Define start and end indices, based on worker count and max points to use
        iSt = (s-1)*nPer+1; iEnd = min(s*nPer, nTime);
        tstru{s} = tstr(iSt:iEnd);
    end
    parfor j = 1:nWorkers
        valcube{j} = ContourTracktest(pre, post, tstr, tstru{j}, tracksFinal, 100, 1, .1, 933, 322, 125,0)
    end
    valcube = cat(2, valcube{:});
    filenameCFP=[pre{1} tstru{1} post{1}];
    save(['data_xy' filenameCFP{1:2} '.mat']);
    displaytext=['-------------------', xystring, ' is done-------------------------'];
    disp(displaytext)
    clear all
end

for h=10:21
    xystring= ['xy', num2str(h)]
    [pre,post,tstr]=getNames(3,1,1,484,1,5,xystring,1,2,'2014-02-15-ekarfra','c1','c2','c3','c4','tif',[10 50]);
    parfor i=1:length(tstr)
        movieInfo{i} = qtetest(pre,post,tstr(i),2,100,14,.6,0);
    end
    movieInfo = cat(2, movieInfo{:});
    scriptTrackGeneral
    %Set the number of workers you want to use
    %You should indicate more workers than you have data elements (I don't think the code is robust to that)
    nWorkers = 9;
    
    %For your application, with tstr (this could be any data that you wanted to splice)
    nTime = numel(tstr);
    
    %Pre-define the cell containing the spliced data
    tstru = cell(nWorkers,1);
    %Calculate the max number of points to put in each cell
    nPer = ceil(nTime/nWorkers);
    
    for s = 1:nWorkers
        %Define start and end indices, based on worker count and max points to use
        iSt = (s-1)*nPer+1; iEnd = min(s*nPer, nTime);
        tstru{s} = tstr(iSt:iEnd);
    end
    parfor j = 1:nWorkers
        valcube{j} = ContourTracktest(pre, post, tstr, tstru{j}, tracksFinal, 100, 1, .1, 933, 322, 125,0)
    end
    valcube = cat(2, valcube{:});
    filenameCFP=[pre{1} tstru{1} post{1}];
    save(['data_xy' filenameCFP{1:2} '.mat']);
    displaytext=['-------------------', xystring, ' is done-------------------------'];
    disp(displaytext)
    clear all
end


