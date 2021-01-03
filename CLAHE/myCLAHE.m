function out = myCLAHE(I, flag)

%ADAPTHISTEQ Contrast-limited Adaptive Histogram Equalization (CLAHE).
%   ADAPTHISTEQ enhances the contrast of images by transforming the
%   values in the intensity image I.  Unlike HISTEQ, it operates on small
%   data regions (tiles), rather than the entire image. Each tile's 
%   contrast is enhanced, so that the histogram of the output region
%   approximately matches the specified histogram. The neighboring tiles 
%   are then combined using bilinear interpolation in order to eliminate
%   artificially induced boundaries.  The contrast, especially
%   in homogeneous areas, can be limited in order to avoid amplifying the
%   noise which might be present in the image.

%--------------------------- The algorithm ----------------------------------
%
%  1. Obtain all the inputs: 
%    * image
%    * number of regions in row and column directions
%    * number of bins for the histograms used in building image transform
%      function (dynamic range)
%    * clip limit for contrast limiting (normalized from 0 to 1)
%    * other miscellaneous options
%  2. Pre-process the inputs:  
%    * determine real clip limit from the normalized value
%    * if necessary, pad the image before splitting it into regions
%  3. Process each contextual region (tile) thus producing gray level mappings
%    * extract a single image region
%    * make a histogram for this region using the specified number of bins
%    * clip the histogram using clip limit
%    * create a mapping (transformation function) for this region
%  4. Interpolate gray level mappings in order to assemble final CLAHE image
%    * extract cluster of four neighboring mapping functions
%    * process image region partly overlapping each of the mapping tiles
%    * extract a single pixel, apply four mappings to that pixel, and 
%      interpolate between the results to obtain the output pixel; repeat
%      over the entire image
%
%  See code for further details.
%
%-----------------------------------------------------------------------------

dimI = size(I);% 图像大小

% 'NumTiles'  Two-element vector of positive integers: [M N].
%                    [M N] specifies the number of tile rows and
%                    columns.  Both M and N must be at least 2. 
%                    The total number of image tiles is equal to M*N.
%                    Default: [8 8].
numTiles = [8 8];
%size of the single tile
dimTile = dimI ./ numTiles;
%Controls the range of the output image data. e.g. [0 255] for uint8
fullRange = [0 255];

%   'NBins'        Positive integer scalar.
%                  Sets number of bins for the histogram used in building a
%                  contrast enhancing transformation. Higher values result 
%                  in greater dynamic range at the cost of slower processing
%                  speed.
%
%                  Default: 256.
numBins = 256;
%   'normClipLimit'    Real scalar from 0 to 1. Default: 0.01.
%   'ClipLimit'             limits contrast enhancement. Higher numbers 
%                                result in more contrast. 
%compute actual clip limit from the normalized value entered by the user
%maximum value of normClipLimit=1 results in standard AHE, i.e. no clipping;
%the minimum value minClipLimit would uniformly distribute the image pixels
%across the entire histogram, which would result in the lowest possible
%contrast value

%prod为乘法计算;ceil为向上取整;round为四舍五入

normClipLimit = 0.01;
numPixInTile = prod(dimTile);
minClipLimit = ceil(numPixInTile/numBins);
clipLimit = minClipLimit + round(normClipLimit*(numPixInTile-minClipLimit));
tileMappings = makeTileMappings(I, numTiles, dimTile, numBins, clipLimit, ...
                                 fullRange);

%Synthesize the output image based on the individual tile mappings. 
if flag
    out = makeClaheImage(I, tileMappings, numTiles, ...
                         dimTile);
else
    out = makeClaheImage2(I, tileMappings, numTiles, ...
                         dimTile);
end
end

%-----------------------------------------------------------------------------

function tileMappings = ...
    makeTileMappings(I, numTiles, dimTile, numBins, clipLimit,...
                      fullRange)%计算得到每一个patch的CLAHE直方图，按对应顺序排列到对应位置.
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

%直方图均衡函数
function mapping = makeMapping(imgHist,  fullRange, ...
    numPixInTile)%得到patch每一个灰度值变换后对应的灰度值，即表现为mapping

histSum = cumsum(imgHist);% cumsum：第1行到第m行的所有元素累加和。histSum即为原图像的累积分布函数
valSpread  = fullRange(2) - fullRange(1);%255-0
  scale =  valSpread/numPixInTile;
  mapping = min(fullRange(1) + histSum*scale,...
                fullRange(2)); %limit to max (列向量)
end

%-----------------------------------------------------------------------------
% This function clips the histogram according to the clipLimit and
% redistributes clipped pixels across bins below the clipLimit
function imgHist = clipHistogram(imgHist, clipLimit, numBins)
% total number of pixels overflowing clip limit in each bin
totalExcess = sum(max(imgHist - clipLimit,0));  %求出在所有灰度值上多于设置的limit阈值的像素点个数

% clip the histogram and redistribute the excess pixels in each bin
avgBinIncr = floor(totalExcess/numBins);%将所有超过门限的像素值相加求平均
upperLimit = clipLimit - avgBinIncr; % bins larger than this will be
                                     % set to clipLimit

% this loop should speed up the operation by putting multiple pixels
% into the "obvious" places first
for k=1:numBins
  if imgHist(k) > clipLimit
    imgHist(k) = clipLimit;
  else
    if imgHist(k) > upperLimit % high bin count
      totalExcess = totalExcess - (clipLimit - imgHist(k));
      imgHist(k) = clipLimit;
    else
      totalExcess = totalExcess - avgBinIncr;
      imgHist(k) = imgHist(k) + avgBinIncr;      
    end
  end
end

% this loops redistributes the remaining pixels, one pixel at a time
k = 1;
while (totalExcess ~= 0)
  %keep increasing the step as fewer and fewer pixels remain for
  %the redistribution (spread them evenly)
  stepSize = max(floor(numBins/totalExcess),1);
  for m=k:stepSize:numBins
    if imgHist(m) < clipLimit
      imgHist(m) = imgHist(m)+1;
      totalExcess = totalExcess - 1; %reduce excess
      if totalExcess == 0
        break;
      end
    end
  end
  if stepSize~=1
    k = k+1; %prevent from always placing the pixels in bin #1
    if k > numBins % start over if numBins was reached
      k = 1;
    end
  end
end

end

%-----------------------------------------------------------------------------
% This function interpolates between neighboring tile mappings to produce a 
% new mapping in order to remove artificially induced tile borders.  
% Otherwise, these borders would become quite visible.  The resulting
% mapping is applied to the input image thus producing a CLAHE processed
% image.

function claheI = makeClaheImage(I, tileMappings, numTiles, ...
                                  dimTile)
%initialize the output image to zeros (preserve the class of the input image)
claheI = I;
claheI(:) = 0;

imgTileRow=1;
for k=1:numTiles(1)+1
  if k == 1  %special case: top row
    imgTileNumRows = dimTile(1)/2; %图像块维数的一半
    mapTileRows = [1 1];
  else 
    if k == numTiles(1)+1 %special case: bottom row      
      imgTileNumRows = dimTile(1)/2;
      mapTileRows = [numTiles(1) numTiles(1)];%【8 8】
    else %default values
      imgTileNumRows = dimTile(1); %图像块维数
      mapTileRows = [k-1, k]; %[upperRow lowerRow]
    end
  end
  
  % loop over columns of the tileMappings cell array
  imgTileCol=1;
  for l=1:numTiles(2)+1
    if l == 1 %special case: left column
      imgTileNumCols = dimTile(2)/2;
      mapTileCols = [1, 1];
    else
      if l == numTiles(2)+1 % special case: right column
        imgTileNumCols = dimTile(2)/2;
        mapTileCols = [numTiles(2), numTiles(2)];
      else %default values
        imgTileNumCols = dimTile(2);
        mapTileCols = [l-1, l]; % right left
      end
    end
    
    % Extract four tile mappings
    ulMapTile = tileMappings{mapTileRows(1), mapTileCols(1)};
    urMapTile = tileMappings{mapTileRows(1), mapTileCols(2)};
    blMapTile = tileMappings{mapTileRows(2), mapTileCols(1)};
    brMapTile = tileMappings{mapTileRows(2), mapTileCols(2)};

    % Calculate the new greylevel assignments of pixels 坐标设置完成
    
    % within a submatrix of the image specified by imgTileIdx. This 
    % is done by a bilinear interpolation between four different mappings 
    % in order to eliminate boundary artifacts.
    
    normFactor = imgTileNumRows*imgTileNumCols; %normalization factor  

    for i = 1:imgTileNumRows
        
        for j = 1: imgTileNumCols
            rowidx = imgTileRow+i-1;
            colidx  = imgTileCol+j-1;
            
            claheI(rowidx, colidx) = ...
                ((imgTileNumRows-i)*((imgTileNumCols-j)*double(ulMapTile(I(rowidx,colidx)+1))+...
                    j*double(urMapTile(I(rowidx,colidx)+1)))+...
                  i*((imgTileNumCols-j)*double(blMapTile(I(rowidx,colidx)+1))+...
                    j*double(brMapTile(I(rowidx,colidx)+1))))...
                 /normFactor;
                
        end
    end
    imgTileCol = imgTileCol + imgTileNumCols;    
  end %over tile cols
  imgTileRow = imgTileRow + imgTileNumRows;
end %over tile rows
end

function claheI = makeClaheImage2(I, tileMappings, numTiles, dimTile)
claheI = I;
claheI(:) = 0;
imgTileRow=1;
for i = 1:numTiles(1)
    imgTileCol=1;
    for j = 1:numTiles(2)
        for row = 1:dimTile(1)     
            for col = 1: dimTile(2)
                rowidx = imgTileRow+row-1;
                colidx = imgTileCol+col-1;
                claheI(rowidx, colidx) = ...
                    tileMappings{i, j}(I(rowidx,colidx)+1);
            end
        end
        imgTileCol = imgTileCol + dimTile(1);
    end
    imgTileRow = imgTileRow + dimTile(2);
end
end
