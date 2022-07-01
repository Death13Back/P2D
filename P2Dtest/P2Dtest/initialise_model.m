function [init_point, n_alg, initial_terminal_voltage] = initialise_model(param)
% initialise_model 执行所有感兴趣的变量的分析初始化



sign_input_density = evaluate_sign_input_density(param); % 根据操作模式评估输入电流/功率密度的符号

T = param.T_init * ones(param.Nsum+param.Nal+param.Ncu,1);

cs_star_p   = param.cs_p_init*ones(param.Np, 1);
cs_star_n   = param.cs_n_init*ones(param.Nn, 1);
cs_star     = [cs_star_p; cs_star_n];

[Phis_p_init,~,Phis_n_init,~] = param.OpenCircuitPotentialFunction(cs_star,T,param,sign_input_density);

Phis_init = [Phis_p_init; Phis_n_init];

if param.OperatingMode==1 || param.OperatingMode==4
    I_density           = param.I_density;
elseif param.OperatingMode==2 || param.OperatingMode==5
    I_density           = param.P_density/(Phis_init(1)-Phis_init(end)); % 从平衡开始不需要线性插值
else
    I_density = 0; % CV 模式的虚拟值
end

ce_init = param.ce_init*[ones(param.Np,1);ones(param.Ns,1);ones(param.Nn,1)];
Phie_init = zeros(param.Np + param.Ns + param.Nn, 1); % 指定为该系统的地电位

solverFlux  = zeros(param.Np+param.Nn, 1);
film = zeros(param.Nn, 1);

jflux_init = ionicFlux(ce_init, cs_star, Phis_init, Phie_init, T, solverFlux, film, param,sign_input_density,I_density);

if param.EnableAgeing == 1
    js_init = 0.483e-5*ones(param.Nn, 1); % 副反应通量的完全任意/随机值
else
    js_init = zeros(param.Nn, 1);
end

init_point = [
    jflux_init;...
    Phis_init;...
    Phie_init;...
    js_init;...
    I_density
    ];

n_alg = length(init_point);
initial_terminal_voltage = diff([Phis_n_init(end);Phis_p_init(1)]); % 从平衡开始不需要线性插值
end
