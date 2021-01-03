function Visualize_SSE(sample, k)
%查看SSE随K值变化而变化的图像
dim = size(sample);
coordinate = zeros(k, 2);
for i = 1:k
    coordinate(i, 1) = i;
    if i == 1
        avg = mean(sample);
        d = zeros(dim(1), 1);
        for j = 1:dim(1)
            X = [sample(j, :);avg];
            d(j,1) = cal_dist(X, 2);
        end
        coordinate(i, 2) = sum(d);
    else
        [~, ~, sumd]=Kmeans(sample, i, 0.1, 9000);
        sumd = sumd.^2;
        coordinate(i, 2)=sum(sumd);
    end
end
figure('name','K-SSE')
plot(coordinate(:, 1),coordinate(:, 2))
xlabel('K')
ylabel('SSE')
end