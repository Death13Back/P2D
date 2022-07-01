function [ce_mean_p, ce_mean_ps, ce_mean_s, ce_mean_sn, ce_mean_n] = interpolateElectrolyteConcentration(ce,param)
%	interpolateElectrolyteConcentration 使用调和平均值在控制体积的边界插入电解质浓度值。



%% 电解质浓度插值

% 正极内插值
beta_ce_p = 0.5;
ce_mean_p = ce(1:param.Np-1).*ce(2:param.Np)./ (beta_ce_p*ce(2:param.Np) + (1-beta_ce_p)*ce(1:param.Np-1));

% 隔膜与正极界面的插值
beta_ce_ps = param.deltax_p*param.len_p/2 / (param.deltax_p*param.len_p/2 + param.deltax_s*param.len_s/2);
ce_mean_ps = ce(param.Np)*ce(param.Np+1)/(beta_ce_ps*ce(param.Np+1) + (1-beta_ce_ps)*ce(param.Np));

% 隔膜内的插值
beta_ce_s = 0.5;
ce_mean_s = ce(param.Np+1:param.Np+param.Ns-1).*ce(param.Np+2:param.Np+param.Ns)./ (beta_ce_s*ce(param.Np+2:param.Np+param.Ns) + (1-beta_ce_s)*ce(param.Np+1:param.Np+param.Ns-1));

% 在隔膜和负极之间的界面上进行插值
beta_ce_sn = param.deltax_s*param.len_s/2 / (param.deltax_n*param.len_n/2 + param.deltax_s*param.len_s/2);
ce_mean_sn = ce(param.Np+param.Ns)*ce(param.Np+param.Ns+1)/(beta_ce_sn*ce(param.Np+param.Ns+1) + (1-beta_ce_sn)*ce(param.Np+param.Ns));

% 负极内插值
beta_ce_n = 0.5;
ce_mean_n = ce(param.Np+param.Ns+1:end-1).*ce(param.Np+param.Ns+2:end)./ (beta_ce_n*ce(param.Np+param.Ns+2:end) + (1-beta_ce_n)*ce(param.Np+param.Ns+1:end-1));

end