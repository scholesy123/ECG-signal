function [time_axis,abnomial,thoracic] = read_data(filename)
%return
%time_axis:2500*1 时间点
%abnomial ：2500*5 5条数据
%thoracic : 2500*3 3条数据
data = importdata(filename); 
time_axis = data(:,1);
abnomial = data(:,2:6);
thoracic = data(:,7:9);
end


