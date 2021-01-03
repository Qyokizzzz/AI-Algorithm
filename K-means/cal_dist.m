function dist = cal_dist(X,p)
%计算两个样本点之间的闵可夫斯基距离
%当p=2时即为欧氏距离
%当p=1时即为曼哈顿距离
%p趋近于无穷时为切比雪夫距离
dist = sum( (X(1,:) - X(2,:)).^2 )^(1/p);
end