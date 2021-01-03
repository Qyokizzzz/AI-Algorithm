function tileMappings = makeTileMappings(I, numTiles, dimTile, numBins, clipLimit, fullRange)
%numTiles:8*8
%dimTile:每一块中的像素数
%numBins:256个色域
%fullRange:[0:255]
numPixInTile = prod(dimTile);%一个patch中的像素个数

tileMappings = cell(numTiles);%创建空单位数组：生成8*8维度的矩阵
% extract and process each tile
imgCol = 1;
for col=1:numTiles(2),
  imgRow = 1;
  for row=1:numTiles(1),
    
    tile = I(imgRow:imgRow+dimTile(1)-1,imgCol:imgCol+dimTile(2)-1);
    % input parsing of imhist
    tileHist = imhist(tile, numBins); % column vector，在这一块区域上的直方图  
    tileHist = clipHistogram(tileHist, clipLimit, numBins);
    tileMapping = makeMapping(tileHist,  fullRange, ...
                              numPixInTile);%直方图均衡映射
    % assemble individual tile mappings by storing them in a cell array;
    tileMappings{row,col} = tileMapping;%结果存入最终的直方图矩阵
    imgRow = imgRow + dimTile(1); %变换到下一个块
  end
  imgCol = imgCol + dimTile(2); % move to the next column of tiles
end

%-----------------------------------------------------------------------------
% Calculate the equalized lookup table (mapping) based on cumulating the input 
% histogram.  
end