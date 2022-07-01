function [res,Up,Un,dudt_p,dudt_n,Keff,J_S] = algebraicStates(x,ce,cs_barrato,Q,T,film,param,t)
% algebraicStates 计算各个单元中心处所有代数方程的残差。



% 离子通量
jflux     = x(param.jflux_indices-param.ndiff);
% 固相电位
Phis      = x(param.Phis_indices-param.ndiff);
% 电解质电位
Phie      = x(param.Phie_indices-param.ndiff);
% 副反应通量
js        = x(param.js_indices-param.ndiff);

I_density = x(param.curr_dens_indices-param.ndiff);

sign_input_density = evaluate_sign_input_density(param); % 根据操作模式评估输入电流/功率密度的符号

if param.OperatingMode==3 % 检查是否在恒电位充电模式下运行
    param.I_density = I_density;
end

% 固体电势上的残差
res_Phis = solidPhasePotential(jflux + [zeros(param.Np,1);js],...
    param,...
    Phis,I_density);

% 电解液电位的残差
[res_Phie, Keff] = electrolytePotential(jflux + [zeros(param.Np,1);js],...
    ce,...
    T,...
    param,...
    Phie);

% 表面平均浓度
cs_star = surfaceConcentration(cs_barrato,...
    jflux,...
    Q,...
    T,...
    param);

% 离子通量计算
[jflux_calc,Up,Un,dudt_p,dudt_n,J_S] = ionicFlux(ce,...
    cs_star,...
    Phis,...
    Phie,...
    T,...
    jflux + [zeros(param.Np,1);js],...
    film,...
    param,sign_input_density,I_density);

%% 构建残差数组
% 离子通量残差
jflux_res = jflux-jflux_calc;
% 副反应残差
js_res    = js-J_S;


switch param.edge_values
	% 通过线性插值
	case 2
		phis_pos_cc = 1.5*Phis(1) - 0.5*Phis(2);
		phis_neg_cc = -0.5*Phis(end-1) + 1.5*Phis(end);
	% 考虑质心值
	otherwise
		phis_pos_cc = Phis(1);
		phis_neg_cc = Phis(end);
end

% 下面的变量“scalar_res”是标量残差变量的占位符，表示电流密度的残差（在模式 1、2、4 和 5 中）或电压残差（在模式 3 中）
if param.OperatingMode==1 || param.OperatingMode==4
    scalar_res = I_density-param.I_density;
    res = [res_Phie;res_Phis;jflux_res;js_res;scalar_res]; % 返回残差（注意 res_Phis 包括根据施加的电流密度的两个电极的边界条件）
elseif param.OperatingMode==2 || param.OperatingMode==5
    scalar_res = I_density-(param.P_density/(phis_pos_cc-phis_neg_cc));
    pwr_dens_neg_BC_res = + (param.P_density) ...
                          + (param.sig_eff(1)/(param.deltax_p*param.len_p))*(Phis(2)*phis_pos_cc-Phis(1)*phis_pos_cc) ...
                          + (param.sig_eff(3)/(param.deltax_n*param.len_n))*(Phis(end-1)*phis_neg_cc-Phis(end)*phis_neg_cc) ...
                          - (param.len_p*param.deltax_p*param.a_i(1)*param.F*jflux(1))*phis_pos_cc ...
                          - (param.len_n*param.deltax_n*param.a_i(3)*param.F*jflux(end))*phis_neg_cc;
    res = [res_Phie;res_Phis;pwr_dens_neg_BC_res;jflux_res;js_res;scalar_res]; % 返回残差（注意负极的功率密度残差已经包含在残差向量中）
else
    scalar_res = phis_pos_cc-phis_neg_cc-param.V_reference;
    res = [res_Phie;res_Phis;jflux_res;js_res;scalar_res]; % 需要验证这一点。 返回残差（对于恒电位情况）
end
end
