function [cs_average_init, ce_init, T_init, film_init, Q_init, n_diff] =   differentialInitialConditions(param)
%	DifferentInitialConditions 初始化微分变量的值



% 检查用于固体扩散的模型类型
if(param.SolidPhaseDiffusion==1 || param.SolidPhaseDiffusion==2)
    % 使用简化模型时使用此初始化
    cs_average_init     = [param.cs_p_init*ones(param.Np,1);param.cs_n_init*ones(param.Nn,1)];
elseif (param.SolidPhaseDiffusion==3)
    % 如果使用完整模型（Fick 定律），则修改初始条件以说明固相扩散方程结构。
    cs_average_init     = [param.cs_p_init*ones(param.Np*param.Nr_p,1);param.cs_n_init*ones(param.Nn*param.Nr_n,1)];
end
% 其他微分变量的初始值。           
ce_init             = param.ce_init*[ones(param.Np,1);ones(param.Ns,1);ones(param.Nn,1)];
T_init              = param.T_init * ones(param.Nsum+param.Nal+param.Ncu,1);
film_init           = zeros(param.Nn,1);
Q_init              = zeros(param.Np+param.Nn,1);

% 在每个cell中存储微分变量的数量。
n_diff= sum([length(cs_average_init) length(ce_init) length(T_init) length(film_init) length(Q_init)]);
end
