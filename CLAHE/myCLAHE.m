function out = myCLAHE(I, flag)
dimI = size(I);
numTiles = [8 8];
dimTile = dimI ./ numTiles;
fullRange = [0 255];
numBins = 256;
normClipLimit = 0.01;
numPixInTile = prod(dimTile);
minClipLimit = ceil(numPixInTile/numBins);
clipLimit = minClipLimit + round(normClipLimit*(numPixInTile-minClipLimit));
tileMappings = makeTileMappings(I, numTiles, dimTile, numBins, clipLimit, ...
                                 fullRange);

if flag
    out = makeClaheImage(I, tileMappings, numTiles, ...
                         dimTile);
else
    out = makeClaheImage2(I, tileMappings, numTiles, ...
                         dimTile);
end
end
