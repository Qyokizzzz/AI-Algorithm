%查看三维聚类效果
load fisheriris
data = normalization(meas); 
X = data(:, 1:3);
[idx, C, ~] = Kmeans(X, 3, 0, Inf);
figure('name', '三维聚类效果')
plot3(C(:, 1), C(:, 2), C(:, 3), 'kx')
hold on
view(3)
plot3(X(idx==1,1),X(idx==1,2),X(idx==1,3),'r*')
hold on
plot3(X(idx==2,1),X(idx==2,2),X(idx==2,3),'b.')
hold on
plot3(X(idx==3,1),X(idx==3,2),X(idx==3,3),'gx')