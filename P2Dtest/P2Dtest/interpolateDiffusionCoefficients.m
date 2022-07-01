function [Deff_p_medio, Deff_s_medio, Deff_n_medio] = interpolateDiffusionCoefficients(Deff_p,Deff_s,Deff_n,param)
% interpolateDiffusionCoefficients 评估扩散系数值的插值。



%% %% 正极的扩散系数
% 与隔膜界面的 Beta 系数
betaD_ps = param.deltax_p*param.len_p/2 / (param.deltax_p*param.len_p/2 + param.deltax_s*param.len_s/2);
%正极内的  Beta 系数
betaD_p = 0.5;
% 正极内系数的调和平均值
Deff_p_medio = (Deff_p(1:end-1).*Deff_p(2:end)./(betaD_p*Deff_p(2:end)+(1-betaD_p)*Deff_p(1:end-1)));
% 与隔膜界面处系数的调和平均值
Deff_medio_ps = Deff_p(end)*Deff_s(1) / (betaD_ps*Deff_s(1)+(1-betaD_ps)*Deff_p(end));
% 构建整个数组
Deff_p_medio = [Deff_p_medio;Deff_medio_ps];
%% 隔膜的扩散系数
% 与负极界面的 Beta 系数
betaD_sn = param.deltax_s*param.len_s/2 / (param.deltax_n*param.len_n/2 + param.deltax_s*param.len_s/2);
%隔膜内的  Beta 系数
betaD_s = 0.5;
% 分离器内系数的 调和平均值
Deff_s_medio = (Deff_s(1:end-1).*Deff_s(2:end))./(betaD_s*Deff_s(2:end)+(1-betaD_s)*Deff_s(1:end-1));
% 界面处系数的调和平均值
% 电极
Deff_medio_sn = Deff_n(1)*Deff_s(end) / (betaD_sn*Deff_n(1)+(1-betaD_sn)*Deff_s(end));
% 构建整个数组
Deff_s_medio = [Deff_s_medio;Deff_medio_sn];
%% 负极的扩散系数
%隔膜内的  Beta 系数
betaD_n = 0.5;
% 构建整个数组。 （注意数组末尾加了一个0
% 以获得正确的尺寸；最后 0 将是
% 替换为 A_tot 矩阵中的正确值）。
Deff_n_medio = [(Deff_n(1:end-1).*Deff_n(2:end))./(betaD_n*Deff_n(2:end)+(1-betaD_n)*Deff_n(1:end-1));0 ];

end
