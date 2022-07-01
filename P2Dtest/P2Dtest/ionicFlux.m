function [jflux,U_p,U_n,dudt_p,dudt_n,J_s] = ionicFlux(ce,cs_star,Phis,Phie,T,solverFlux,film,param,sign_input_density,I_density)
% onicFlux 计算锂离子在电极-电解质界面处的摩尔通量密度 [mol/(m^2*s)]。



%% 正极
% 计算正极和负极的 OCV。
[U_p,dudt_p,U_n,dudt_n] = param.OpenCircuitPotentialFunction(cs_star,T,param,sign_input_density);

% 计算反应速率。
[k_pT, k_nT] = param.ReactionRatesFunction(T,param);

% 正电极离子通量
deltap = ((0.5*param.F)./(param.R*T(param.Nal+1:param.Nal+param.Np))).*(Phis(1:param.Np)-Phie(1:param.Np)-U_p);
ip = 2*k_pT.*sqrt(ce(1:param.Np)).*sqrt(cs_star(1:param.Np)).*sqrt(param.cs_maxp-cs_star(1:param.Np));
jnp_calc = ip.* sinh(deltap);

%% 负极

% 如果启用老化，请考虑 SEI 电阻
if(param.EnableAgeing==1)
    eta_n   = (Phis(param.Np+1:end)-Phie(param.Np+param.Ns+1:end)-U_n -param.F*solverFlux(param.Np+1:end).*(param.R_SEI+film./(param.k_n_aging)));
else
    eta_n   = (Phis(param.Np+1:end)-Phie(param.Np+param.Ns+1:end)-U_n);
end

deltan      = ((0.5*param.F)./(param.R*T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn))).*eta_n;
in          = 2*k_nT.*sqrt(ce(param.Np+param.Ns+1:end)).*sqrt(cs_star(param.Np+1:end)).*sqrt(param.cs_maxn-cs_star(param.Np+1:end));
jnn_calc    = in.* sinh(deltan);

J_s = zeros(param.Nn,1);

% % 当施加的电流密度是符号变量时切换情况
% if(isa(I_density,'casadi.MX')||isa(I_density,'casadi.SX') && param.EnableAgeing == 1)
%     eta_s = Phis(param.Np+1:end) - Phie(param.Np+param.Ns+1:end) - param.Uref_s - param.F*solverFlux(param.Np+1:end).*(param.R_SEI+film./(param.k_n_aging));
%     % 副反应通量的塔菲尔方程。
%     alpha   = 0.5*param.F./(param.R*T(param.Nal+param.Np+param.Ns+1:end-param.Ncu));
%     % 借助CasADi的if_else语句，可以表示根据符号量I_密度的值进行切换的动力学
%     J_s = if_else(I_density>=0,-param.i_0_jside.*(I_density/param.I1C)^param.w.*(exp(-alpha.*eta_s))./param.F,zeros(param.Nn,1));
% elseif(param.EnableAgeing == 1 && sign_input_density > 0)
%     eta_s = Phis(param.Np+1:end) - Phie(param.Np+param.Ns+1:end) - param.Uref_s - param.F*solverFlux(param.Np+1:end).*(param.R_SEI+film./(param.k_n_aging));
%     % 副反应通量的塔菲尔方程。
%     alpha   = 0.5*param.F./(param.R*T(param.Nal+param.Np+param.Ns+1:end-param.Ncu));
%     J_s = -param.i_0_jside.*(I_density/param.I1C)^param.w.*(exp(-alpha.*eta_s))./param.F;
% end
%% 返回值
jflux = [jnp_calc;jnn_calc];
