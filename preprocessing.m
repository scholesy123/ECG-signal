function[new_data]=preprocessing(data)
%250Hz 10s
%高通3hz截止
%低通巴特沃斯65hz截止
%陷波中心频率50hz
data = highpass(data,3,250);
[y,x] = butter(2,65/125);
data = filter(y,x,data);
w=50/(250/2);
bw=w;
[num,den]=iirnotch(w,bw); % notch filter implementation 
new_data=filter(num,den,data);
%new_data = rescale(new_data,-1,1);
end