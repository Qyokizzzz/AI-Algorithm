function claheI = makeClaheImage(I, tileMappings, numTiles,dimTile)
%插值
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