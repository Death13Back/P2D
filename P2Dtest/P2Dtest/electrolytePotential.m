function [res_Phie,Keff] = electrolytePotential(jflux,ce,T,param,Phie)
% 电解质电位评估空间离散的电解质电位方程的残差（线法）。



%% 有效电解质电导率
% 出于基准测试目的对此进行修改
Keff_p = param.ElectrolyteConductivityFunction(ce(1:param.Np),T(param.Nal+1:param.Nal+param.Np),param,'p');
Keff_s = param.ElectrolyteConductivityFunction(ce(param.Np+1:param.Np+param.Ns),T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns),param,'s');
Keff_n = param.ElectrolyteConductivityFunction(ce(param.Np+param.Ns+1:end),T(param.Nal+param.Np+param.Ns+1:end-param.Ncu),param,'n');

% 出于基准测试目的取消注释
% Keff_p = param.eps_p^param.brugg_p *(4.1253*1e-2 + 5.007*1e-4*ce(1:param.Np) - 4.7212*1e-7*ce(1:param.Np).^2 +1.5094*1e-10*ce(1:param.Np).^3 -1.6018*1e-14*ce(1:param.Np).^4);
% Keff_s = param.eps_s^param.brugg_s *(4.1253*1e-2 + 5.007*1e-4*ce(param.Np+1:param.Np+param.Ns) - 4.7212*1e-7*ce(param.Np+1:param.Np+param.Ns).^2 +1.5094*1e-10*ce(param.Np+1:param.Np+param.Ns).^3 -1.6018*1e-14*ce(param.Np+1:param.Np+param.Ns).^4);
% Keff_n = param.eps_n^param.brugg_n *(4.1253*1e-2 + 5.007*1e-4*ce(param.Np+param.Ns+1:end) - 4.7212*1e-7*ce(param.Np+param.Ns+1:end).^2 +1.5094*1e-10*ce(param.Np+param.Ns+1:end).^3 -1.6018*1e-14*ce(param.Np+param.Ns+1:end).^4);

Keff = [Keff_p;Keff_s;Keff_n];

% 由于 Keff 的值是在每个 CV 的中心计算的，因此需要对这些量进行插值并在 CV 的边缘找到它们的值
[Keff_p_medio, Keff_s_medio, Keff_n_medio] = interpolateElectrolyteConductivities(Keff_p,Keff_s,Keff_n,param);

%% 矩阵构建

% 第 i 个元素
A_p = diag(Keff_p_medio);
A_p(2:end,2:end)    = A_p(2:end,2:end) + diag(Keff_p_medio(1:end-1));
% 第 i+1 个元素
A_p(1:end-1,2:end)  = A_p(1:end-1,2:end) - diag(Keff_p_medio(1:end-1));
% 第 i-1 个元素
A_p(2:end,1:end-1)  = A_p(2:end,1:end-1) - diag(Keff_p_medio(1:end-1));

% 第 i 个元素
A_s = diag(Keff_s_medio);
A_s(2:end,2:end)    = A_s(2:end,2:end) + diag(Keff_s_medio(1:end-1));
% 第 i+1 个元素
A_s(1:end-1,2:end)  = A_s(1:end-1,2:end) - diag(Keff_s_medio(1:end-1));
% 第 i-1 个元素
A_s(2:end,1:end-1)  = A_s(2:end,1:end-1) - diag(Keff_s_medio(1:end-1));

% 第 i 个元素
A_n = diag(Keff_n_medio);
A_n(2:end,2:end)    = A_n(2:end,2:end) + diag(Keff_n_medio(1:end-1));
% 第 i+1 个元素
A_n(1:end-1,2:end)  = A_n(1:end-1,2:end) - diag(Keff_n_medio(1:end-1));
% 第 i-1 个元素
A_n(2:end,1:end-1)  = A_n(2:end,1:end-1) - diag(Keff_n_medio(1:end-1));

A_n = A_n./(param.deltax_n*param.len_n);
A_s = A_s./(param.deltax_s*param.len_s);
A_p = A_p./(param.deltax_p*param.len_p);

A_tot = blockDiagonalMatrix(param,A_p,A_s,A_n);

% 固定值以在正极左侧强制执行 BC。
% A_tot(1,1:2) = [Keff_p_medio(1) -Keff_p_medio(1)]./(param.deltax_p*param.len_p);

% 负极最后体积的 Phie 值是已知且固定的。
% 现在，我们有 -k_eff（最后一个内部面）+ keff（最后一个内部面）
% 只保留下面两行之一的注释，以便在负极末端强制执行 phi_e 边界条件

switch param.edge_values
	case 2
		% 在负电极 = 0 处设置 phi_e 的间接方式。（线性插值）
		A_tot(end,end-1:end) = [-1/3 1];  
	otherwise
		% 在负电极处设置 phi_e 的直接方式 = 0（整个 CV 设置为 0）
		A_tot(end,end-1:end) = [0 1]; 
end

%% Interfaces 正极（正极的最后体积）

% 这里我们在最后一个体积元
den_sn = (param.deltax_p*param.len_p/2+param.deltax_s*param.len_s/2);
last_p = Keff_p_medio(end-1)/(param.deltax_p*param.len_p);
A_tot(param.Np,param.Np-1:param.Np+1) = [-last_p (last_p+Keff_p_medio(end)/den_sn) -Keff_p_medio(end)/den_sn];
%% Interfaces 正极（隔膜的第一体积）

% 这里我们在隔膜的第一体积单元
den_sp 	= (param.deltax_p*param.len_p/2+param.deltax_s*param.len_s/2);
first_s = Keff_s_medio(1)/(param.deltax_s*param.len_s);
A_tot(param.Np+1,param.Np:param.Np+2) = [-Keff_p_medio(end)/den_sp (first_s+Keff_p_medio(end)/den_sp) -first_s];

%% Interfaces 正极界面（隔膜的最后体积单元）
% 这里我们在隔膜的最后一个体积单元
den_sn = (param.deltax_n*param.len_n/2+param.deltax_s*param.len_s/2);
last_s = Keff_s_medio(end-1)/(param.deltax_s*param.len_s);
A_tot(param.Np+param.Ns,param.Np+param.Ns-1:param.Np+param.Ns+1) = [-last_s (last_s+Keff_s_medio(end)/den_sn) -Keff_s_medio(end)/den_sn];

%% Interfaces 正极界面（负极的第一体积单元）
% 这里我们在负极第一体积单元里面
den_ns 	= (param.deltax_n*param.len_n/2+param.deltax_s*param.len_s/2);
first_n = Keff_n_medio(1)/(param.deltax_n*param.len_n);
A_tot(param.Np+param.Ns+1,param.Np+param.Ns:param.Np+param.Ns+2) = [-Keff_s_medio(end)/den_ns (first_n+Keff_s_medio(end)/den_sn) -first_n];

% 电解质浓度插值 评估控制体积边缘处电解质浓度值的插值。
[ce_mean_p, ce_mean_ps, ce_mean_s, ce_mean_sn, ce_mean_n] = interpolateElectrolyteConcentration(ce,param);
%% 温度插值
% 评估控制体积边缘的温度值
[T_mean_p, T_mean_ps, T_mean_s, T_mean_sn, T_mean_n] = interpolateTemperature(T,param);
%% 电解质通量
% 评估电解质浓度通量的插值
% 控制体积的边缘。
[ce_flux_p, ce_flux_ps, ce_flux_s, ce_flux_sn, ce_flux_n] = interpolateElectrolyteConcetrationFluxes(ce,param);
%% RHS 阵列
K = 2*param.R*(1-param.tplus) / param.F;

prod_p = (Keff_p_medio.*[T_mean_p;T_mean_ps].*[ce_flux_p;ce_flux_ps].*[1./ce_mean_p;1/ce_mean_ps]);

prod_s = (Keff_s_medio.*[T_mean_s;T_mean_sn].*[ce_flux_s;ce_flux_sn].*[1./ce_mean_s;1/ce_mean_sn]);

% 不需要 Keff_n_medio 的最后一个单元
prod_n = (Keff_n_medio(1:end-1).*T_mean_n.*ce_flux_n.*1./ce_mean_n);

if(isa(prod_p,'casadi.SX') || isa(prod_p,'casadi.MX'))
    prod_tot = [prod_p;prod_s;prod_n];
    prod_tot = prod_tot(2:end)-prod_tot(1:end-1);
else
    prod_tot = diff([prod_p;prod_s;prod_n]);
end
prod_tot = [prod_p(1);prod_tot];

flux_p = param.deltax_p*param.len_p*param.F*param.a_i(1)*jflux(1:param.Np);
flux_s = param.deltax_s*param.len_s*param.F*param.a_i(2);
flux_n = param.deltax_n*param.len_n*param.F*param.a_i(3)*jflux(param.Np+1:end-1);

flux_tot = [flux_p;flux_s*ones(param.Ns,1);flux_n];

f = flux_tot-K*prod_tot;

% 将 Phie 的最后一个元素设置为 0（强制执行 BC）
f=[f;0];

if(~isa(prod_p,'casadi.SX') && ~isa(prod_p,'casadi.MX'))
    A_tot = sparse(A_tot);
end

% 返回电解液电位的残值
res_Phie = A_tot*Phie-f;
end
