function [U_p,dudt_p,U_n,dudt_n] = openCircuitPotential(cs_star,T,param,inputCurrentSign)
% openCircuitPotential 评估特定电池阴极和阳极的开路电压。
% 计量单位为 [V]
%
% [U_p,dudt_p,U_n,dudt_n] = openCircuitPotential(cs_star,T,param,inputCurrentSign)
%
%       - 输出：
% - U_p 和 U_n ：分别代表阳极和阴极的 OCV。
% - dudt_p 和 dudt_n ：分别表示相对于温度变化的 OCV 变化。
% 这仅在启用温度的模拟中考虑。
%
%
% 您可以修改 OCV 的计算方式，作为电解质浓度和温度的函数。 主脚本还将传递 param 数组。
	



% 计算电池正极的开路电压
	theta_p  = cs_star(1:param.Np)./param.cs_maxp;

	% 计算 OCV 相对于温度变化的变化 [V/K]
	dudt_p   = -0.001 * (0.199521039-0.928373822*theta_p + 1.364550689000003*theta_p.^2-0.6115448939999998*theta_p.^3);
	dudt_p   = dudt_p./(1-5.661479886999997*theta_p +11.47636191*theta_p.^2-9.82431213599998*theta_p.^3+3.048755063*theta_p.^4);

	% 定义正极的 OCV
	U_p    = (-4.656+88.669*theta_p.^2 - 401.119*theta_p.^4 + 342.909*theta_p.^6 - 462.471*theta_p.^8 + 433.434*theta_p.^10);
	U_p      = U_p./(-1+18.933*theta_p.^2-79.532*theta_p.^4+37.311*theta_p.^6-73.083*theta_p.^8+95.96*theta_p.^10);
	U_p      = U_p + (T(param.Nal+1:param.Nal+param.Np)-param.Tref).*dudt_p;



% 计算电池在负极的开路电压 
	theta_n  = cs_star(param.Np+1:end)./ param.cs_maxn;

	% 计算 OCV 相对于温度变化的变化 [V/K]
	dudt_n = 0.001*(0.005269056 +3.299265709*theta_n-91.79325798*theta_n.^2+1004.911008*theta_n.^3-5812.278127*theta_n.^4 + ...
		19329.7549*theta_n.^5 - 37147.8947*theta_n.^6 + 38379.18127*theta_n.^7-16515.05308*theta_n.^8); % There is a typo in this eqn. in the original paper
	dudt_n = dudt_n./(1-48.09287227*theta_n+1017.234804*theta_n.^2-10481.80419*theta_n.^3+59431.3*theta_n.^4-195881.6488*theta_n.^5+...
		374577.3152*theta_n.^6 - 385821.1607*theta_n.^7 + 165705.8597*theta_n.^8);

	% 定义负极的 OCV
	U_n   = 0.7222 + 0.1387*theta_n + 0.029*theta_n.^0.5 - 0.0172./theta_n + 0.0019./theta_n.^1.5 + 0.2808*exp(0.9-15*theta_n)-0.7984*exp(0.4465*theta_n - 0.4108);
	U_n   = U_n + (T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)-param.Tref).*dudt_n;
end