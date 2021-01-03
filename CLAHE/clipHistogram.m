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