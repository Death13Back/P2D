function [Qrev, Qrxn, Qohm]=heatGenerationRates(Phis,Phie,jflux,T,cs_star,ce,param)
% heatGenerationRates 评估热模型中使用的热源项。



sign_input_density = evaluate_sign_input_density(param); % 根据操作模式评估输入电流/功率密度的符号

% 检索有效的电解质电导系数。
Keff_p = param.ElectrolyteConductivityFunction(ce(1:param.Np),T(param.Nal+1:param.Nal+param.Np),param,'p');
Keff_s = param.ElectrolyteConductivityFunction(ce(param.Np+1:param.Np+param.Ns),T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns),param,'s');
Keff_n = param.ElectrolyteConductivityFunction(ce(param.Np+param.Ns+1:end),T(param.Nal+param.Np+param.Ns+1:end-param.Ncu),param,'n');

% 评估 Qohm 计算中使用的导数
[dPhis, dPhie, dCe] = ThermalDerivatives(Phis',Phie',ce',param);	% y(end) 从包含代数状态的向量获得电流
dPhis = dPhis';
dPhie = dPhie';
dCe = dCe';
[Up,dudt_p,Un,dudt_n] = param.OpenCircuitPotentialFunction(cs_star,T,param,sign_input_density);

%%可逆发热率

% 正极
Qrev_p = param.F*param.a_i(1)*jflux(1:param.Np).*T(param.Nal+1:param.Nal+param.Np).*dudt_p;

% 负极
Qrev_n = param.F*param.a_i(3)*jflux(param.Np+1:end).*T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn).*dudt_n;

%% 反应生热率

% 正过电位
eta_p = (Phis(1:param.Np)-Phie(1:param.Np)-Up);
% 正反应放热率
Qrxn_p = param.F*param.a_i(1)*jflux(1:param.Np).*eta_p;

% 负过电位
eta_n = (Phis(param.Np+1:end)-Phie(param.Np+param.Ns+1:end)-Un);
% 负反应生热率
Qrxn_n = param.F*param.a_i(3)*jflux(param.Np+1:end).*eta_n;

%% 欧姆发热率

% 正极欧姆发生率
Qohm_p = param.sig_eff(1) * (dPhis(1:param.Np)).^2 + Keff_p.*(dPhie(1:param.Np)).^2 + 2*param.R*Keff_p.*T(param.Nal+1:param.Nal+param.Np)*(1-param.tplus)/param.F.*dCe(1:param.Np).*1./ce(1:param.Np).*dPhie(1:param.Np);
% 隔膜欧姆发生率
Qohm_s = Keff_s.*(dPhie(param.Np+1:param.Np+param.Ns)).^2 + 2*param.R*Keff_s.*T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns)*(1-param.tplus)/param.F.*dCe(param.Np+1:param.Np+param.Ns).*1./ce(param.Np+1:param.Np+param.Ns).*dPhie(param.Np+1:param.Np+param.Ns);
% 负电极欧姆发生率
Qohm_n = param.sig_eff(3) * (dPhis(param.Np+1:end)).^2 +Keff_n.*(dPhie(param.Np+param.Ns+1:end)).^2 + 2*param.R*Keff_n.*T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)*(1-param.tplus)/param.F.*dCe(param.Np+param.Ns+1:end).*1./ce(param.Np+param.Ns+1:end).*dPhie(param.Np+param.Ns+1:end);

Qrev = [Qrev_p zeros(1,param.Ns) Qrev_n];
Qrxn = [Qrxn_p zeros(1,param.Ns) Qrxn_n];
Qohm = [Qohm_p Qohm_s Qohm_n];
