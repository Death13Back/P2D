function [Dps_eff, Dns_eff] = solidPhaseDiffusionCoefficients(T,param)
% solidPhaseDiffusionCoefficients 计算固相的扩散系数 [m^2 /s]。
% 用户可以修改脚本以满足特定要求。
Dps_eff     = param.Dps*ones(param.Np,1);
Dns_eff     = param.Dns*ones(param.Nn,1);

% if(param.TemperatureEnabled>=1)
%     Dps_eff     = param.Dps*exp(-param.EaDps/param.R*(1./T(param.Nal+1:param.Nal+param.Np)-1/param.Tref));
% else
%     Dps_eff     = param.Dps*ones(param.Np,1);
% end
% 
% if(param.TemperatureEnabled>=1)
%     Dns_eff     = param.Dns*exp(-param.EaDns/param.R*(1./T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)-1/param.Tref));
% else
%     Dns_eff     = param.Dns*ones(param.Nn,1);
% end
% 
% end
