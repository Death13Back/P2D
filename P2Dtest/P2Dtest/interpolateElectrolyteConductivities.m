function [Keff_p_medio, Keff_s_medio, Keff_n_medio] = interpolateElectrolyteConductivities(Keff_p,Keff_s,Keff_n,param)
%	interpolateElectrolyteConductivitys 使用调和平均值在控制体积的边界插入电解质电导率。



%% 正极平均电导率
beta_p = 0.5;
Keff_p_medio = Keff_p(1:end-1).*Keff_p(2:end)./(beta_p*Keff_p(2:end)+(1-beta_p)*Keff_p(1:end-1));

% Keff_p_medio 是在正极-隔膜的调和平均值

beta_p_s = param.deltax_p*param.len_p/2 /(param.deltax_s*param.len_s/2+param.deltax_p*param.len_p/2);

Keff_p_s_interface = Keff_p(end)*Keff_s(1) / (beta_p_s*Keff_s(1) + (1-beta_p_s)*Keff_p(end));

Keff_p_medio = [Keff_p_medio;Keff_p_s_interface];

%% 隔膜平均电导率
%计算隔膜的调和平均值
beta_s = 0.5;
Keff_s_medio = Keff_s(1:end-1).*Keff_s(2:end)./(beta_s*Keff_s(2:end)+(1-beta_s)*Keff_s(1:end-1));

% Keff_s_medio 是在隔膜与负极的调和平均值

beta_s_n = param.deltax_s*param.len_s/2 /(param.deltax_s*param.len_s/2+param.deltax_n*param.len_n/2);

Keff_s_n_interface = Keff_s(end)*Keff_n(1) / (beta_s_n*Keff_n(1) + (1-beta_s_n)*Keff_s(end));

Keff_s_medio = [Keff_s_medio;Keff_s_n_interface];

%% 负极平均电导率
% 计算负极的调和平均值

beta_n = 0.5;
Keff_n_medio = Keff_n(1:end-1).*Keff_n(2:end)./(beta_n*Keff_n(2:end)+(1-beta_n)*Keff_n(1:end-1));

Keff_n_medio = [Keff_n_medio;0]; % 不使用cc 接口。 零只是为了匹配尺寸
end