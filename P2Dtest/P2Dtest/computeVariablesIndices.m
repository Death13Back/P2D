function param = computeVariablesIndices(param)
%   computeVariablesIndices 计算代码中使用的微分变量的索引



param.ce_indices         = (1:param.Nsum);

% 根据使用的模型修改固相指数。 如果使用菲克定律，则还需要考虑固体颗粒内部的扩散。
if(param.SolidPhaseDiffusion==1 || param.SolidPhaseDiffusion==2)
    param.cs_average_indices    = (param.ce_indices(end)+1:param.ce_indices(end)+param.Np+param.Nn);
elseif(param.SolidPhaseDiffusion==3)
    param.cs_average_indices = (param.ce_indices(end)+1:param.ce_indices(end)+param.Np*param.Nr_p+param.Nn*param.Nr_n);
end

param.T_indices    = (param.cs_average_indices(end)+1:param.cs_average_indices(end)+param.Nal+param.Nsum+param.Ncu);
param.film_indices = (param.T_indices(end)+1:param.T_indices(end)+param.Nn);
param.Q_indices    = (param.film_indices(end)+1:param.film_indices(end)+param.Np+param.Nn);

% 存储代数变量的索引。
param.jflux_indices     = (param.Q_indices(end)+1:param.Q_indices(end)+param.Np+param.Nn);
param.Phis_indices      = (param.jflux_indices(end)+1:param.jflux_indices(end)+param.Np+param.Nn);
param.Phie_indices      = (param.Phis_indices(end)+1:param.Phis_indices(end)+param.Np+param.Ns+param.Nn);
param.js_indices        = (param.Phie_indices(end)+1:param.Phie_indices(end)+param.Nn);
param.curr_dens_indices = (param.js_indices(end)+1);

end
