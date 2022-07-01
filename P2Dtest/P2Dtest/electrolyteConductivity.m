function Keff = electrolyteConductivity(ce,T,param,batterySection)
%electrodeConductivity 评估电解质相的电导系数。 计量单位为[S/m]
%
% Keff =electrodeConductivity(ce,T,param) 评估电池的阳极、隔膜和阴极电解质相的电导系数。 您可以修改脚本以满足您的特定需求。
%
% 可以在等温情况（param.TemperatureEnabled=0）或绝热情况（param.TemperatureEnabled=1 或 2）下评估电导系数。
%
% 您可以修改电导系数的计算方式，作为电解质浓度和温度的函数。 主脚本还将传递 param 数组。
        switch(batterySection)
        case'p'
            Keff = param.eps_p^param.brugg_p *(4.1253*1e-2 + 5.007*1e-4*ce - 4.7212*1e-7*ce.^2 +1.5094*1e-10*ce.^3 -1.6018*1e-14*ce.^4);
        case 's'
            Keff = param.eps_s^param.brugg_s *(4.1253*1e-2 + 5.007*1e-4*ce - 4.7212*1e-7*ce.^2 +1.5094*1e-10*ce.^3 -1.6018*1e-14*ce.^4);
        case 'n'
            Keff = param.eps_n^param.brugg_n *(4.1253*1e-2 + 5.007*1e-4*ce - 4.7212*1e-7*ce.^2 +1.5094*1e-10*ce.^3 -1.6018*1e-14*ce.^4);
        end


% if(param.TemperatureEnabled>=1)
%     switch(batterySection)
%         case'p'
%             Keff = param.eps_p^param.brugg_p *(1e-4*ce.*((-10.5+0.668*1e-3*ce+0.494*1e-6*ce.^2) +...
%                 (0.074  -1.78*1e-5*ce -8.86*1e-10*ce.^2).*T + (-6.96*1e-5+2.8*1e-8*ce).*T.^2).^2);
%         case 's'
%             Keff = param.eps_s^param.brugg_s *(1e-4*ce.*((-10.5+0.668*1e-3*ce+0.494*1e-6*ce.^2) +...
%                 (0.074  -1.78*1e-5*ce -8.86*1e-10*ce.^2).*T + (-6.96*1e-5+2.8*1e-8*ce).*T.^2).^2);
%         case 'n'
%             Keff = param.eps_n^param.brugg_n *(1e-4*ce.*((-10.5+0.668*1e-3*ce+0.494*1e-6*ce.^2) +...
%                 (0.074  -1.78*1e-5*ce -8.86*1e-10*ce.^2).*T + (-6.96*1e-5+2.8*1e-8*ce).*T.^2).^2);
%     end
% else
%         switch(batterySection)
%         case'p'
%             Keff = param.eps_p^param.brugg_p *(4.1253*1e-2 + 5.007*1e-4*ce - 4.7212*1e-7*ce.^2 +1.5094*1e-10*ce.^3 -1.6018*1e-14*ce.^4);
%         case 's'
%             Keff = param.eps_s^param.brugg_s *(4.1253*1e-2 + 5.007*1e-4*ce - 4.7212*1e-7*ce.^2 +1.5094*1e-10*ce.^3 -1.6018*1e-14*ce.^4);
%         case 'n'
%             Keff = param.eps_n^param.brugg_n *(4.1253*1e-2 + 5.007*1e-4*ce - 4.7212*1e-7*ce.^2 +1.5094*1e-10*ce.^3 -1.6018*1e-14*ce.^4);
%         end
% end