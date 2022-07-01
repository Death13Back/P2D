function cs_star = surfaceConcentration(cs_barrato,jflux,Q,T,param)
% surfaceConcentration 评估电极表面锂离子的浓度。



% 固相扩散系数
[Dps_eff, Dns_eff] = param.SolidDiffusionCoefficientsFunction(T,param);

% 检查选择了哪种固体扩散模型。
if(param.SolidPhaseDiffusion==1) % 二参数模型
    % 评估两个电极中的平均表面浓度。
    % 阴极
    cs_star_p = cs_barrato(1:param.Np)-(param.Rp_p./(Dps_eff.*5)).*jflux(1:param.Np);
    % 阳极
    cs_star_n = cs_barrato(param.Np+1:end)-(param.Rp_n./(Dns_eff.*5)).*jflux(param.Np+1:end);
elseif(param.SolidPhaseDiffusion==2) % 三参数模型
    % 阴极
    cs_star_p = cs_barrato(1:param.Np)+(param.Rp_p./(Dps_eff.*35)).*(-jflux(1:param.Np)+8*Dps_eff.*Q(1:param.Np));
    % 阳极
    cs_star_n = cs_barrato(param.Np+1:end)+(param.Rp_n./(Dns_eff.*35)).*(-jflux(param.Np+1:end)+8*Dns_eff.*Q(param.Np+1:end));
elseif(param.SolidPhaseDiffusion==3) % 完整模型
    p_indices = param.Nr_p:param.Nr_p:param.Nr_p*param.Np;
    n_indices = param.Nr_n:param.Nr_n:param.Nr_n*param.Nn;
    % 如果使用完整模型，只需取r=Rp处的浓度数据
    cs_star_p = cs_barrato(p_indices);
    cs_star_n = cs_barrato(n_indices+p_indices(end));
end
% 返回残差
cs_star = [cs_star_p;cs_star_n];

end
