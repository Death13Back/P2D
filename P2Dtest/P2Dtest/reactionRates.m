function [k_pT, k_nT] = reactionRates(T,param)
% REACTIONRATES 阴极和阳极的反应速率 (k) [m^2.5/(m^0.5 s)]。
% 用户可以修改此脚本以满足特定要求。
k_pT     = param.k_p;
k_nT     = param.k_n;

% if(param.TemperatureEnabled>=1)
%     k_pT     = param.k_p*exp(-param.Eakip/param.R*(1./T(param.Nal+1:param.Nal+param.Np)-1/param.Tref));
% else
%     k_pT     = param.k_p;
% end
% 
% if(param.TemperatureEnabled>=1)
%     k_nT     = param.k_n*exp(-param.Eakin/param.R*(1./T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)-1/param.Tref));
% else
%     k_nT     = param.k_n;
% end
% 
% end
