function mapping = makeMapping(imgHist,fullRange,numPixInTile)
%得到映射关系
histSum = cumsum(imgHist);
%cumsum：第1行到第m行的所有元素累加和。histSum即为原图像的累积分布函数
valSpread  = fullRange(2) - fullRange(1);%255-0
  scale =  valSpread/numPixInTile;
  mapping = min(fullRange(1)+histSum*scale,fullRange(2));
end