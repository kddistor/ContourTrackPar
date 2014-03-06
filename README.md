CountourTracking
================

Robust cell tracking. Nuclear and Cyto region intensity grabber.

####Usage:

[pre,post,tstr]=getNames(3,1,1,194,1,2,'xy09',1,2,'2012-09-25-ekar2','c1','c2','c3','c4','tif',[16 50]);

movieInfo=qte(pre,post,tstr,2,[16 50]);

scriptTrackGeneral

[valcube cbound nbound] = ContourTrack(pre, post, tstr, tracksFinal, [7 30])

overlayTracksObjectsMovieC(pre, post, tstr, 1, cbound, nbound, tracksFinal, 1:534)


figure, imagesc(valcube(:,:,7)) to display nuclear ratios

figure, imagesc(valcube(:,:,7)) to display cyto ratios
