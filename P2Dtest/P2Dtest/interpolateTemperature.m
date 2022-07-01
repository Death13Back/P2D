function [T_mean_p, T_mean_ps, T_mean_s, T_mean_sn, T_mean_n] = interpolateTemperature(T,param)
%	interpolateTemperature 使用调和平均值评估控制体积边缘的温度插值。



% 正极内插补
beta_T_p = 0.5;
T_mean_p = T(param.Nal+1:param.Nal+param.Np-1).*T(param.Nal+2:param.Nal+param.Np)./ (beta_T_p*T(param.Nal+2:param.Nal+param.Np) + (1-beta_T_p)*T(param.Nal+1:param.Nal+param.Np-1));

% 隔膜与正极界面的插值
beta_T_ps = param.deltax_p*param.len_p/2 / (param.deltax_p*param.len_p/2 + param.deltax_s*param.len_s/2);
T_mean_ps = T(param.Nal+param.Np)*T(param.Nal+param.Np+1)/(beta_T_ps*T(param.Nal+param.Np+1) + (1-beta_T_ps)*T(param.Nal+param.Np));

% 隔膜内的插值
beta_T_s = 0.5;
T_mean_s = T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns-1).*T(param.Nal+param.Np+2:param.Nal+param.Np+param.Ns)./ (beta_T_s*T(param.Nal+param.Np+2:param.Nal+param.Np+param.Ns) + (1-beta_T_s)*T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns-1));

% 在隔膜和负极之间的界面上进行插值
beta_T_sn = param.deltax_s*param.len_s/2 / (param.deltax_n*param.len_n/2 + param.deltax_s*param.len_s/2);
T_mean_sn = T(param.Nal+param.Np+param.Ns)*T(param.Nal+param.Np+param.Ns+1)/(beta_T_sn*T(param.Nal+param.Np+param.Ns+1) + (1-beta_T_sn)*T(param.Nal+param.Np+param.Ns));

% 负极内插值
beta_T_n = 0.5;
T_mean_n = T(param.Nal+param.Np+param.Ns+1:end-(param.Ncu)-1).*T(param.Nal+param.Np+param.Ns+2:end-param.Ncu)./ (beta_T_n*T(param.Nal+param.Np+param.Ns+2:end-param.Ncu) + (1-beta_T_n)*T(param.Nal+param.Np+param.Ns+1:end-(param.Ncu+1)));

end