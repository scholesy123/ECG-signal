
% [w, y, e] = CLMS(2, 0.0000015, 0.00000017, 400, 4, w1_init, w2_init, 0, x, d);
filename = "FOETAL_ECG.txt";
[time_axis,abnomial,thoracic] = read_data(filename);
d = abnomial(:,2);
x = thoracic(:,3);
mu = 0.0017;
M=6;
delta = 1e-7;
x_new = preprocessing(x);
d_new = preprocessing(d);
[yn, W, en1]=nlmsFunc(x_new, d_new, M, mu, delta);
%[en2,w,y] =RLSFilterIt(x_new,d_new);
rlsFilt = dsp.RLSFilter(6,'ForgettingFactor', 0.95);
[y,en2] = rlsFilt(x_new,d_new);
%set(gcf,'position',[3 3 8 2]);%figture位置，最下角，宽高
%set (gca,'position',[0.1,0.5,0.9,0.3] );%axis位置，最下角，宽高
%plot(x_new);
%figure;plot(d_new);
order = 6;
w1_init = ones(order, 1)*0.005;
w2_init = ones(order, 1)*0.003;
[w3, y3, e3, all_w3] = CLMS(2, 1e-7, 1e-8, 400, 4, w1_init, w2_init, 0, x_new, d_new);
[w2, y2, e2, all_w2] = CRLS(40000, 400, 4, 0.995, 0.95, w1_init, w2_init, 0, x_new, d_new);
[w1, y1, e1, all_w1] = LMSPLUSRLS(40, 1e-7,0.995, 400, 4, w1_init, w2_init, 0, x_new, d_new);
%SNR
SNRe1 = 10*log10(sum(e1.^2)/sum(y1.^2));
SNRe2 = 10*log10(sum(e2.^2)/sum(y2.^2));
SNRe3 = 10*log10(sum(e3.^2)/sum(y3.^2));
SNR = 10*log10(sum(en2.^2)/sum(y.^2));
SNRen = 10*log10(sum(en1.^2)/sum(yn.^2));
%[w22, y22, e22, all_w22] = CRLS(1, 400, 4, 1, 1, w1_init, w2_init, 0, y2, d_new);
%权重收敛
%figure;
%set(gcf,'position',[3 3 8 2]);
%plot(all_w2(:,1:15:2494)');
%plot(all_w2');
%xlabel('Samples');
%ylabel('Weights of filter');
%title('CRLS weights convergence');
% 
% figure;
% %set(gcf,'position',[3 3 8 2]);
% plot(all_w2(:,1:15:2494)');
% %plot(all_w2');
% xlabel('Samples');
% ylabel('Weights of filter');
% title('CRLS weights convergence');
% %波形图

% set(gcf,'position',[3 3 8 2]);
% plot(e22(1600:2000));
% hold on;
% title('res of CLMS');
% %EMSE曲线

figure;
plot(d_new(1600:2400));hold on;plot(e1(1600:2400));
set(gca, 'linewidth', 1.1, 'fontsize', 25, 'fontname', 'times')
xlabel('Samples');
ylabel('Amplitude');
%title('AECG and FECG extracted by RLS-LMS');
legend('AECG','FECG');
figure;
plot(d_new(1600:2400));hold on;plot(e3(1600:2400));
set(gca, 'linewidth', 1.1, 'fontsize', 25, 'fontname', 'times')
xlabel('Samples');
ylabel('Amplitude');
%title('AECG and FECG extracted by CLMS');
legend('AECG','FECG');
figure;
plot(d_new(1600:2400));hold on;plot(en2(1600:2400));
set(gca, 'linewidth', 1.1, 'fontsize', 25, 'fontname', 'times')
xlabel('Samples');
ylabel('Amplitude');
%title('AECG and FECG extracted by RLS');
legend('AECG','FECG');
figure;
plot(d_new(1600:2400));hold on;plot(en1(1600:2400));
set(gca, 'linewidth', 1.1, 'fontsize', 25, 'fontname', 'times')
xlabel('Samples');
ylabel('Amplitude');
%title('AECG and FECG extracted by NLMS');
legend('AECG','FECG');
figure;
plot(d_new(1600:2400))
set(gca, 'linewidth', 1.1, 'fontsize', 25, 'fontname', 'times')
xlabel('Samples');
ylabel('Amplitude');
%title('AECG samples from 1600 to 2400');
figure;
plot(d_new(1600:2400));hold on;plot(e2(1600:2400));
set(gca, 'linewidth', 1.1, 'fontsize', 25, 'fontname', 'times')
xlabel('Samples');
ylabel('Amplitude');
%title('AECG and FECG extracted by CRLS');
legend('AECG','FECG');



%subplot(2,3,2),plot(d_new(1600:2400));
%xlabel('samples');
%ylabel('Amplitude');
%title('AECG');
%subplot(2,3,3),plot(d_new(1600:2400));hold on;plot(e1(1600:2400));
%subplot(2,3,4),plot(d_new(1600:2400));hold on;plot(e3(1600:2400));
%subplot(2,3,5),plot(d_new(1600:2400));hold on;plot(en1(1600:2400));
%subplot(2,3,6),plot(d_new(1600:2400));hold on;plot(en2(1600:2400));






                
