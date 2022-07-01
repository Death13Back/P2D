function [res_Phis] = solidPhasePotential(jflux,param,Phis,I_density)
% solidPhasePotential 计算固体电位方程的残差。



%% 正极

% RHS 为正极中的固体电位。 强制执行左侧的 BC (Neumann BC)
f_p = ((param.len_p*param.deltax_p*param.a_i(1)*param.F*jflux(1))-I_density)*param.deltax_p*param.len_p/param.sig_eff(1);
% RHS 为正极中的固体电位。
f_p = [f_p;(param.len_p^2*param.deltax_p^2*param.a_i(1)*param.F*jflux(2:param.Np))/param.sig_eff(1)];

%% 负极

% RHS 表示负极中的固体电位。
f_n = (param.len_n^2*param.deltax_n^2*param.a_i(3)*param.F*jflux(param.Np+1:end-1))/param.sig_eff(3);

if param.OperatingMode==1 || param.OperatingMode==4 || param.OperatingMode==3 %右侧的 Neumann BC 仅在使用外加电流密度作为输入进行操作时强制执行
    % RHS 表示负极中的固体电位。
    f_n =[f_n;((param.len_n*param.deltax_n*param.a_i(3)*param.F*jflux(end))+I_density)*param.deltax_n*param.len_n/param.sig_eff(3)];
end

%% 残差数组
% 返回残差数组
res_Phis = [
    param.A_p*Phis(1:param.Np)-f_p;...
    param.A_n*Phis(param.Np+1:end)-f_n;
    ];
end
