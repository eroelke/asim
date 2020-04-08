function table = reformat_data(data,len)

table = nan(len,size(data,2));

x_min = data(1,1);
x_max = data(end,1);
dx = (x_max-x_min)/(len-1);
x = x_min:dx:x_max;

for ii = 1:length(x)
    table(ii,1) = x(ii);
    table(ii,2:end) = interp1(data(:,1), data(:,2:end), x(ii), 'linear');
end

return