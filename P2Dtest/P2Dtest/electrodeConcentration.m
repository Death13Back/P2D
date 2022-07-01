function [ddCs, rhsCs] = electrodeConcentration(dCs,cs_barrato,T,jflux,param)
%电极浓度描述了电极内锂离子浓度的 ODE。
if(param.SolidPhaseDiffusion==1 || param.SolidPhaseDiffusion==2)
    % 阴极
    rhsCs_p  =((-3/param.Rp_p)*jflux(1:param.Np));
    ddCs_p   = dCs(1:param.Np) - rhsCs_p;
    
    % 阳极
    rhsCs_n  = ((-3/param.Rp_n)*jflux(param.Np+1:end));
    ddCs_n   = dCs(param.Np+1:end) - rhsCs_n;
else
    switch param.SolidPhaseDiffusionNumericalScheme
        % 使用 FDM 方法进行固相扩散
        case 1
            [rhsCs_p, rhsCs_n, ddCs_p, ddCs_n] = FDM9orderElectrodeDiffusion(T, cs_barrato, jflux, dCs, param);
        % 使用谱方法离散化固相扩散
        case 2
            [rhsCs_p, rhsCs_n, ddCs_p, ddCs_n] = spectralMethodElectrodeDiffusion(T, cs_barrato, jflux, dCs, param);
    end
end
rhsCs   = [rhsCs_p;rhsCs_n];
ddCs    = [ddCs_p;ddCs_n];
end
