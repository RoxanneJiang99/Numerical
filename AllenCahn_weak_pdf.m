% % pdf of tamed Allen Cahn
clc; 
clear; clear all; 

randn('state',100);
% 参数设定
rounds = 5000; % number of paths sampled
N = 2^5; %spatial discretization




T = 100; %final time point
tau = 2^(-4); % temporal stepsize

%T = 5; %final time point
%tau = 2^(-6); % temporal stepsize


M = T / tau; % temporal discretization

%beta1=6;   %Example1   Version1:13
beta2=1;   %Example2  Version1:22

A = -pi^2*(1:N).^2 ;  
kappa = 0; %噪声类型：kappa=0为I-Wiener;kappa=2为Q-Wiener;
 
%Wiener Process:
Q = ones(1,N); 
for i =1 : N 
    Q(i) = 1 / (1 + i.*(log(i)).^2);
end



YInitial1 = ones(1,N);
YInitial2 = zeros(1,N);
Y0 = zeros(1,N);
for i =1:N
       Y0(i) = exp(sqrt(N));
end
YInitial3 = Y0 ;
Y0 = zeros(1,N);
for i =1:N
        Y0(i)=3.*pi./2;
end
YInitial4 = Y0 ;
Y0 = zeros(1,N);
for i =1:N
       Y0(i)=3*pi./2 + log(i);
end
YInitial5 = Y0 ;

Y1_pdf = ones( rounds,1 ); %建立一个矩阵，存储每个round最终步的测试函数下的值,用来绘制分布函数
Y2_pdf = ones( rounds,1 ); %建立一个矩阵，存储每个round最终步的测试函数下的值,用来绘制分布函数
Y3_pdf = ones( rounds,1 ); %建立一个矩阵，存储每个round最终步的测试函数下的值,用来绘制分布函数
Y4_pdf = ones( rounds,1 ); %建立一个矩阵，存储每个round最终步的测试函数下的值,用来绘制分布函数
Y5_pdf = ones( rounds,1 ); %建立一个矩阵，存储每个round最终步的测试函数下的值,用来绘制分布函数

for i  = 1 : rounds % 路径循环开始  
     tic % 计时开始  
     Rd = randn(M,N);
     Y = YInitial1; % 比较值的初始值分量
    
     for m = 1 : M        
             y = dst(Y) * sqrt(2);
             %fy = 10.* y - y.^3;
             %fy = ( 10.* y -  y.^3)./( 1 + beta1.* tau.* y.^2 ); %example1
             fy = ( 10.* y -  y.^3)./( ( 1 + (tau + (1./(pi^2.* N^2))  ).* y.^8 ).^(1/4)  );%example2
             fY = idst(fy) / sqrt(2);
             Y = exp(A * tau) .* Y...
               + tau * exp(A * tau) .* fY ...
               + sqrt(Q) .* exp(A * tau) * sqrt(tau) .* Rd(m,:);
     end 
   Y1_pdf(i,1) = Y(1); 
  
   
        Rd = randn(M,N);
     Y = YInitial2; % 比较值的初始值分量
    
     for m = 1 : M        
             y = dst(Y) * sqrt(2);
             %fy = 10.* y - y.^3;
             %fy = ( 10.* y -  y.^3)./( 1 + beta1.* tau.* y.^2 ); %example1
             fy = ( 10.* y -  y.^3)./( ( 1 + (tau + (1./(pi^2.* N^2))  ).* y.^8 ).^(1/4)  );%example2
             fY = idst(fy) / sqrt(2);
             Y = exp(A * tau) .* Y...
               + tau * exp(A * tau) .* fY ...
               + sqrt(Q) .* exp(A * tau) * sqrt(tau) .* Rd(m,:);
     end 
   Y2_pdf(i,1) = Y(1); 
   
   
           Rd = randn(M,N);
     Y = YInitial3; % 比较值的初始值分量
    
     for m = 1 : M        
             y = dst(Y) * sqrt(2);
             %fy = 10.* y - y.^3;
             %fy = ( 10.* y -  y.^3)./( 1 + beta1.* tau.* y.^2 ); %example1
             fy = ( 10.* y -  y.^3)./( ( 1 + (tau + (1./(pi^2.* N^2))  ).* y.^8 ).^(1/4)  );%example2
             fY = idst(fy) / sqrt(2);
             Y = exp(A * tau) .* Y...
               + tau * exp(A * tau) .* fY ...
               + sqrt(Q) .* exp(A * tau) * sqrt(tau) .* Rd(m,:);
     end 
   Y3_pdf(i,1) = Y(1); 
   
   
           Rd = randn(M,N);
     Y = YInitial4; % 比较值的初始值分量
    
     for m = 1 : M        
             y = dst(Y) * sqrt(2);
             %fy = 10.* y - y.^3;
             %fy = ( 10.* y -  y.^3)./( 1 + beta1.* tau.* y.^2 ); %example1
             fy = ( 10.* y -  y.^3)./( ( 1 + (tau + (1./(pi^2.* N^2))  ).* y.^8 ).^(1/4)  );%example2
             fY = idst(fy) / sqrt(2);
             Y = exp(A * tau) .* Y...
               + tau * exp(A * tau) .* fY ...
               + sqrt(Q) .* exp(A * tau) * sqrt(tau) .* Rd(m,:);
     end 
   Y4_pdf(i,1) = Y(1); 
   
   
           Rd = randn(M,N);
     Y = YInitial5; % 比较值的初始值分量
    
     for m = 1 : M        
             y = dst(Y) * sqrt(2);
             %fy = 10.* y - y.^3;
             %fy = ( 10.* y -  y.^3)./( 1 + beta1.* tau.* y.^2 ); %example1
             fy = ( 10.* y -  y.^3)./( ( 1 + (tau + (1./(pi^2.* N^2))  ).* y.^8 ).^(1/4)  );%example2
             fY = idst(fy) / sqrt(2);
             Y = exp(A * tau) .* Y...
               + tau * exp(A * tau) .* fY ...
               + sqrt(Q) .* exp(A * tau) * sqrt(tau) .* Rd(m,:);
     end 
   Y5_pdf(i,1) = Y(1); 
   
   
   
 if mod(i,10) == 0
   disp(['Round: ', num2str(i)])
   toc % 计时终止

end 
end% 路径循环终止

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
    
[y1,x1] = ksdensity( Y1_pdf );
[y2,x2] = ksdensity( Y2_pdf );
[y3,x3] = ksdensity( Y3_pdf );
[y4,x4] = ksdensity( Y4_pdf );
[y5,x5] = ksdensity( Y5_pdf );



plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,'LineWidth',1);
hold on;
set(gca,'FontSize',13);
axis tight
legend({'Initial 1','Initial 2','Initial 3','Initial 4','Initial 5'},'Location','NorthEast');
title('Sampling for the space-time full-discretization schemes at T=100','FontSize',13);

  
