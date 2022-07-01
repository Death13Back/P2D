function [ce_t, cs_barrato_t, T_t, jflux_t, Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t, Q_t,t_tot] = retrieveData(ce_t, cs_barrato_t, T_t, jflux_t, Phis_t, Phie_t,cs_star_t, SOC_t, film_t, js_t, Q_t,t_tot, y, t, param)
%	retrieveData 检索结果结构中返回的数据。



%在积分过程后提取微分变量
ce_t            = [ce_t;y(param.ce_indices)];
cs_barrato_t    = [cs_barrato_t;y(param.cs_average_indices)];
T_t             = [T_t;y(param.T_indices)];
film_t          = [film_t;y(param.film_indices)];
Q_t             = [Q_t;y(param.Q_indices)];
% 积分过程后提取代数变量
jflux_t         = [jflux_t;y(param.jflux_indices)];
Phis_t          = [Phis_t;y(param.Phis_indices)];
Phie_t          = [Phie_t;y(param.Phie_indices)];
js_t            = [js_t;y(param.js_indices)];

% 检查是否使用了菲克扩散定律。 这是定义评估 SOC 的正确方法所必需的。
if(param.SolidPhaseDiffusion~=3)
    cs_average = sum(cs_barrato_t(end,param.Np+1:end))/(param.Nn);  % 负电极中的 cs_average
else
    cs_average = sum(cs_barrato_t(end, (param.Np*param.Nr_p) +1:end))/(param.Nn*param.Nr_n); % 负电极中的 cs_average
end

Sout = 100*((cs_average/param.cs_maxn) - param.theta_min_neg)/(param.theta_max_neg- param.theta_min_neg); % 电池-SOC百分比

cs_star_t       = [cs_star_t;surfaceConcentration(cs_barrato_t(end,:)',jflux_t(end,:)',Q_t(end,:)',T_t(end,:)',param)'];
SOC_t           = [SOC_t;Sout];
t_tot           = [t_tot;t];
end
