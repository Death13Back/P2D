
clear all
clc
%addpath(genpath('I:\P2DSolver\P2D_2\P2Dtest'))



% 定义积分时间。
t0 = 0;
tf = 10^4;
% 定义参数结构。
param{1} = Parameters_init;

% 禁用热动力学
%param{1}.TemperatureEnabled = 0;

% 开始模拟。 请注意，最终积分时间为 10^4，当达到 2.5V 的截止电压时，程序将自动停止。

out1 = startSimulation(t0,tf,[],-15,param);

% 存储雅可比矩阵以供将来计算。 这是可能的，因为不同的场景共享相同的模型结构，唯一不同的是施加的电流密度。
param{1}.JacobianFunction = out1.JacobianFun;

% 运行模拟
out2 = startSimulation(t0,tf,[],-30,param);

% 运行模拟
out3 = startSimulation(t0,tf,[],-60,param);

%% 绘图

figure(1)
plot(out1.time{1},out1.Voltage{1},'LineWidth',1)
hold on
plot(out2.time{1},out2.Voltage{1},'--','LineWidth',1)
plot(out3.time{1},out3.Voltage{1},'-.','LineWidth',1)
xlabel('Time [s]')
ylabel('Voltage [V]')
grid on
box on
