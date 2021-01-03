function [index,C,sumd] = Kmeans(sample, k, threshold, n)
%K均值算法
%C:k个簇中心
%index:聚类后每个样本的标记
%sumd:样本点到相应的簇心的距离
%sample:需要进行聚类的样本
%k:划分簇的个数
%threshold:差异度阈值
%n最大迭代次数
iter = 0;
dim = size(sample);
index = zeros(dim(1), 1);
dist = zeros(k, 1);
C = sample(randperm(dim(1), k), :);
while 1
    sumd = zeros(dim(1), 1);
    for i = 1:dim(1)
        for j = 1:k
            X = [sample(i, :);C(j, :)];
            dist(j) = cal_dist(X, 2);
        end
        [d, idx] = min(dist);
        sumd(i) = d;
        index(i) = idx;
    end
    new_C = zeros(k, dim(2));
    c = 0;
    for i = 1:k
        count = 0;
        for j = 1:dim(1)
            if index(j) == i
                count = count + 1;
                new_C(i, :) = new_C(i, :) + sample(j, :);
            end
        end
        new_C(i, :) = new_C(i, :) / count;
        Y = [new_C(i, :);C(i, :)];
        if cal_dist(Y, 2) <= threshold
            c = c + 1;
        end
    end
    iter = iter + 1;
    if c == k
        break;
    elseif iter > n
        break;
    else
        C = new_C;
    end
end
end