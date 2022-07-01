function Deff = electrolyteDiffusionCoefficients(ce,T,param,batterySection)
%electrodeDiffusionCoefficients 计算电解质相的扩散系数 [m^2/s]。 
%[Deff_p, Deff_s, Deff_n] =electrodeDiffusionCoefficients(ce,T,param) 
%评估电池的阳极、隔膜和阴极的扩散系数。
%
% 注意这是主程序的接口。 作者建议保留函数的名称及其签名，同时只修改脚本的主体。
%
% 可以在等温情况（param.TemperatureEnabled=0）或绝热情况（param.TemperatureEnabled=1 或 2）下评估扩散系数。
%
% 您可以修改计算扩散系数的方式，作为电解质浓度和温度的函数。 主脚本还将传递 param 数组。
%
% 用户可以修改此脚本以满足特定要求。
    switch(batterySection)
        case 'p'
            Deff = repmat(param.Dp*param.eps_p^param.brugg_p,param.Np,1);
        case 's'
            Deff = repmat(param.Ds*param.eps_s^param.brugg_s,param.Ns,1);
        case 'n'
            Deff = repmat(param.Dn*param.eps_n^param.brugg_n,param.Nn,1);
    end

% if(param.TemperatureEnabled>=1)
%     switch(batterySection)
%         case 'p'
%             Deff = param.eps_p^param.brugg_p*1e-4*10.^((-4.43-54./(T-229-5e-3*ce)-0.22e-3*ce));
%         case 's'
%             Deff = param.eps_s^param.brugg_s*1e-4*10.^((-4.43-54./(T-229-5e-3*ce)-0.22e-3*ce));
%         case 'n'
%             Deff = param.eps_n^param.brugg_n*1e-4*10.^((-4.43-54./(T-229-5e-3*ce)-0.22e-3*ce));
%     end
% else
%     switch(batterySection)
%         case 'p'
%             Deff = repmat(param.Dp*param.eps_p^param.brugg_p,param.Np,1);
%         case 's'
%             Deff = repmat(param.Ds*param.eps_s^param.brugg_s,param.Ns,1);
%         case 'n'
%             Deff = repmat(param.Dn*param.eps_n^param.brugg_n,param.Nn,1);
%     end
% end
% 
% end
