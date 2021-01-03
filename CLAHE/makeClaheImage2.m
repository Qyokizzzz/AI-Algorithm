function claheI = makeClaheImage2(I, tileMappings, numTiles, dimTile)
%≤ª≤Â÷µ
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