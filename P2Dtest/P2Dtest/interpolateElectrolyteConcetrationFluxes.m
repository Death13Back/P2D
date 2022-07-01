function [ce_flux_p, ce_flux_ps, ce_flux_s, ce_flux_sn, ce_flux_n] = interpolateElectrolyteConcetrationFluxes(ce,param)
%	interpolateElectrolyteConcetrationFluxes 使用调和平均值在控制体积的边界插入电解质浓度通量。



% 正极内的通量
ce_flux_p = (ce(2:param.Np)-ce(1:param.Np-1))/(param.deltax_p*param.len_p);

% 隔膜正极界面处的通量
ce_flux_ps = (ce(param.Np+1)-ce(param.Np)) / ((param.deltax_p*param.len_p/2+param.deltax_s*param.len_s/2));

% 隔膜内通量
ce_flux_s = (ce(param.Np+2:param.Np+param.Ns)-ce(param.Np+1:param.Np+param.Ns-1))/(param.deltax_s*param.len_s);

% 隔膜-负极界面处的通量
ce_flux_sn = (ce(param.Np+param.Ns+1)-ce(param.Np+param.Ns)) / ((param.deltax_n*param.len_n/2+param.deltax_s*param.len_s/2));

% 负极内的通量
ce_flux_n = (ce(param.Np+param.Ns+2:end)-ce(param.Np+param.Ns+1:end-1))/(param.deltax_n*param.len_n);

end