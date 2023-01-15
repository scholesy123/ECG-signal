function [yn, W, en]=nlmsFunc(xn, dn, M, mu, delta)
itr = length(xn);
en = zeros(itr,1);            
W  = zeros(M,itr);    % 每一列代表-次迭代,初始为0
% 迭代计算
for k = M:itr                  % 第k次迭代
    x = xn(k:-1:k-M+1);        % 滤波器M个抽头的输入
    y = W(:,k-1).' * x;        % 滤波器的输出
    en(k) = dn(k)-y;        % 第k次迭代的误差
    % 滤波器权值计算的迭代式
    W(:,k) = W(:,k-1) + mu*en(k)*x/(delta+x.' * x);
    %W(:,k) = W(:,k-1) + 2*mu*en(k)*x;
end
% for k=M:itr
%     x=xn(k:-1:k-M+1);%滤波器M个抽头输入值
%     yn(k)=w(:,k).'*x;%滤波器输出值
%     e(k)=dn(k)-yn(k);%估计误差
%     w(:,k+1)=w(:,k)+2*u*e(k)*x;%滤波器权系数更新
%     mse(k)=sum(abs(e).^2)/k;%均方误差
% end

%yn = inf * ones(size(xn)); % 初值为无穷大是绘图使用，无穷大处不会绘图
for k = M:itr
    x = xn(k:-1:k-M+1);
    yn(k) = W(:,end).'* x;  % 最终输出结果

end
end