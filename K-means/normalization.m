function data = normalization(x)
%归一化处理数据
d = size(x);
data = zeros(d(1),d(2));
for i = 1:d(2)
    column = x(:,i);
    data(:,i) = (column - min(column))...
        /(max(column)-min(column)); 
end
end