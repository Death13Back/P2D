function param= Parameters_init(init_cell_soc_percent)
%	Parameters_init 定义仿真中使用的参数。 要更改单元/模拟器参数，请在此处修改值.
% 如果未提供百分比形式的初始 SOC，则使用默认值 85.51%
if(nargin<1)
	init_cell_soc_percent = 85.51; % 模拟开始时的 Cell-SOC [%] 以与现有 LIONSIMBA 代码兼容 (<=1.023)
end

% 模拟操作模式

param.OperatingMode = 1;

% 有效选项是：
% 1 - 初始时间和最终时间之间的恒定输入电流密度（默认模式）。
% 2 - 初始时间和最终时间之间的恒定输入功率密度。
% 3 - 恒电位电荷。 在这种情况下，电流被认为是一个变量，充电是在恒定电位下进行的。
% 4 - 在 getInputCurrent 脚本中描述的可变电流密度曲线作为 t 的函数
% 5 - 在 getPowerCurrent 脚本中描述的可变功率密度曲线作为 t 的函数

% 启用或禁用温度动态

% param.TemperatureEnabled = 1;

% 有效选项是：
% 0 - 不考虑热动力学。 运行等温模拟
% 1 - 用完整的 PDE 模型模拟的热动力学
% 2 - 使用降阶集总模型模拟的热动力学


% 如果 param.TemperatureEnabled=2，选择（集总）热模型
% param.lumped_thermal_version = 1;

% 有效选项是：
% 1 - 集总热模型，仅从 eta*current 产生热量
% 2 - 由熵和 eta*current 产生热量的集总热模型

% 从“startSimulation”抑制命令窗口输出的选项
param.suppress_status_prints = 0;   %  1 = 抑制输出。 0 = 正常输出

% 化学计量限制
param.theta_max_pos = 0.49550;  % at 100% cell SOC
param.theta_max_neg = 0.85510;  % at 100% cell SOC
param.theta_min_pos = 0.99174;  % at 0% cell SOC
param.theta_min_neg = 0.01429;  % at 0% cell SOC

%% 通用常数
% 法拉第常数  [C/mol]
param.F         = 96487;
% 气体常数      [J / (mol K)]
param.R         = 8.314;

%% 电池外部几何形状（电池包）
param.pouch_length 	= 332.74e-3;         		% [m] （长尺寸）电池包的长度
param.pouch_width  	= 99.06e-3;          		% [m] (短尺寸) 电池包宽度
tab_width 			= 40e-3;       				% 目前用于电流收集和冷却（表面冷却将导致更小的标签）
tab_length 			= 0.75*param.pouch_width;   % 约 70%-80% 的电池包尺寸，其中长度 (L) 代表较短的电池包尺寸
param.tab_area 		= 2*tab_width*tab_length;   % 有 2 个可用于对流散热的选项卡（仅用于选项卡冷却方法）

%% 切片厚度 [m]
% 铝集流体
param.len_al= 10e-6;
% 正极
param.len_p = 80e-6;
% 隔膜
param.len_s = 25e-6;
% 负极
param.len_n = 88e-6;
% 铜集流体
param.len_cu= 10e-6;
% 电池包
param.len_pouch = 160e-6; % 参考: "Li-Ion Pouch Cells for Vehicle Applications—Studies of Water Transmission and Packing Materials", Pontus Svens, Maria Hellqvist Kjell, Carl Tengstedt, Göran Flodberg and Göran Lindbergh, Energies 2013, 6, 400-410; doi:10.3390/en6010400

%%电池电气/内部几何/设计方面（与标准文献进行比较）
% param.i_1C_density_Northrop_cell 下面是耗尽（Northrop）电池开始所需的放电电流密度
% 在 100%C 下至 LCO 电池的（任意）截止电压为 2.7V，具有此处描述的参数（Northrop电池）。
% 这代表一个电化学层（Al-Pos-Neg-Sep-Cu 组合）的容量的间接测量（由典型文献中发表的原始 Newman 模型模拟。）
param.i_1C_density = 29.23; % [A/m^2]
param.no_of_layers_Northrop_cell = 49;     % 假设没有软件 参考电池（Northrop 电池）的层数，通过手工计算可装入 10 毫米电池的层数计算获得。
param.t_stack = param.no_of_layers_Northrop_cell*(param.len_p + param.len_s + param.len_n) + (ceil(0.5*(param.no_of_layers_Northrop_cell + 1))*param.len_cu) + ceil(0.5*param.no_of_layers_Northrop_cell)*param.len_al; % 可用于填充单元电池（即堆叠层）的长度/厚度
assumed_cell_capacity_Ah = 60; % 假设 Ah 由多个层构成的单元格（每个层都具有此文件中给出的参数）。

% 注意：代码不模拟上述容量的电池。 事实上，软件仅模拟一个电化学层。
% 上述变量仅用于计算所有电化学层的总表面积（用于集总热参数，在用户级脚本中可能有用）
I_1C_cell_amps = assumed_cell_capacity_Ah; % 根据 C-rate 的定义。 （即 60Ah 电池的 1C 电流为 60A）

% 下面计算的变量表示在标准 Newman 模型的典型 1D 离散化中垂直于厚度方向的平面区域）
param.overall_surface_area_for_given_layers = I_1C_cell_amps/param.i_1C_density;  % [m^2] 所有层提供的总表面积

%% 热导率 [ W / (m K) ]

% 铝集流体
param.Lambda_al = 237;  % 来自材料数据表和科学标准
% 正极
param.Lambda_p  = 2.1;
% 隔膜
param.Lambda_s  = 0.16;
% 负极
param.Lambda_n  = 1.7;
% 铜集流体
param.Lambda_cu = 401;  % 来自材料数据表和科学标准

%% 电解质扩散系数 [m^2/s]

% 正极区域
param.Dp = 7.5e-10;
% 隔膜
param.Ds = 7.5e-10;
% 负极区域
param.Dn = 7.5e-10;

%% 密度 [kg / m^3 ]
% 铝集流体
param.rho_al = 2700; % 在室温下。 来自材料数据表和科学标准
% 正极
param.rho_p  = 2500;
% 隔膜
param.rho_s  = 1100; % 我们并不确切知道所使用的隔板材料。 常见的隔膜材料的密度为1070 - 1200 kg/m^3（1070、1100、1200等）。 使用中值
% 负极
param.rho_n  = 2500;
% 铜集流体
param.rho_cu = 8940; % 在室温下。 来自材料数据表和科学标准

% LiPF6电解液
param.rho_LiPF6 = 1290; %  来自 'Characterization of Lithium-Ion Battery Thermal Abuse Behavior Using Experimental and Computational Analysis',doi: 10.1149/2.0751510jes, J. Electrochem. Soc. 2015 volume 162, issue 10, A2163-A2173,表 1 也可以 'Thermal analysis of a cylindrical lithium-ion battery', Xiongwen Zhang, Electrochimica Acta,  56 (2011) 1246–1255, 表 3
% 电池包装材料
param.rho_pouch = 1150; % 'Modeling for the scale-up of a lithium-ion polymer battery',Ui Seong Kim, Chee Burm Shin, Chi-Su Kim, Journal of Power Sources, 2008
% 填料/粘合剂
param.rho_pvdf = 1750;  %来自 表1 'Characterization of Lithium-Ion Battery Thermal Abuse Behavior using Experimental and Computational Analysis',doi: 10.1149/2.0751510jes, J. Electrochem. Soc. 2015 volume 162, issue 10, A2163-A2173

%% 温度设置
param.Tref   = 25 + 273.15; % 环境（环境）温度 [K]

% 参考: 'Reciprocating air flow for Li-ion battery thermal management to improve temperature uniformity" , Rajib Mahamud, Chanwoo Park, Journal of Power Sources, 196 (2011) 5685–5696
param.Tmax   = 55 + 273.15; % 绝对最大允许温度 [K] 温度上限，运行期间任何电池（电池组中）的任何点

%% 比热容 [ J / (kg K) ]
param.Cpal   = 897; % 铝集流体
param.Cpp    = 700; % 正极
param.Cps    = 700; % 隔膜
param.Cpn    = 700; % 负极
param.Cpcu   = 385; % 铜集流体
param.CpLiPF6 =  134.1; % 电解液
% 假设：忽略粘合剂/填料的 Cp，因为它们的含量可以忽略不计
% 此外，外包装在 Cp 计算中也被忽略（但在质量计算中被考虑在内）
param.Cppouch = 1464.8; % 基于小包装材料成分的加权计算


%% 集流体电导率 [S/m]
param.sig_al = 3.55e7;
param.sig_cu = 5.96e7;

%% 电解质孔隙率指数（即电解质体积分数）
% 正极
param.eps_p     = 0.385;
% 隔膜
param.eps_s     = 0.724;
% 负极
param.eps_n     = 0.485;


%% 体积分数
param.eps_fi    = [0.025;0;0.0326];

%% Bruggeman 布鲁格曼系数

% 正极
param.brugg_p   = 4;
% 隔膜
param.brugg_s   = 4;
% 负极
param.brugg_n   = 4;

%% 固相扩散系数 [m^2 / s]

% 正极
param.Dps       = 1e-14;
% 负极
param.Dns       = 3.9e-14;

%% 颗粒表面积 [m^2 / m^3]
% 正极
a_p       = 885000;
% 隔膜
a_s       = 0;
%负极
a_n       = 723600;
% 不要删除。 它在代码中使用
param.a_i       = [a_p;a_s;a_n];

%% Transference 转移数 - 不适用于正/负电极
param.tplus     = 0.364;

%% 反应速率常数 [ m ^ 2.5 / (mol^0.5 s ) ]
% 正极
param.k_p       = 2.334e-11;
% 隔膜
param.k_s       = 0;
% 负极
param.k_n       = 5.031e-11;

%% 热交换系数 [W/m^2 K]
param.hcell     = 1;                    % 仅用于沿厚度方向的一维热模型

%% 固相中锂离子的最大浓度 [ mol/m^3 ]
% 正极
param.cs_maxp   = 51554;
% 隔膜
param.cs_maxs   = 0;
% 负极
param.cs_maxn   = 30555;

% 固体颗粒半径 [m] - 正负极相等
param.Rp_p     = 2e-6;
param.Rp_n     = 2e-6;

% 固相电导率 (S/m)
param.sig    = [100;... % 正极
    0;...   % 隔膜
    100     % 负极
    ];

param.vol_fraction_solidphase = (1 - [param.eps_p;param.eps_s;param.eps_n] - param.eps_fi);
param.vol_fraction_solidphase(2) = 0;   % 隔膜中无固相材料
% 有效固相电导率 (S/m)
param.sig_eff = param.sig.*param.vol_fraction_solidphase;

%% 集总热参数
% [param,surface_area_per_face_for_49_layers,~,~] = compute_lumped_mass_and_Cp_avg_for_given_layer_fcn(param.no_of_layers_Northrop_cell,param); % 计算单元的质量和 Cp_avg 并将它们附加到参数结构（用于集总热模型计算）
% param.surface_area_per_face_Northrop_cell = surface_area_per_face_for_49_layers;   % 这将是映射到 10 毫米厚的特定 WXL 电池的 49 层诺斯罗普Northrop电池所特有的每个面的表面积This will be the surface area per face specific to the 49 layer Northrop cell mapped to a 10mm thick pouch of spefic WXL
% param.h_lumped = 150;    % [W/(m^2 K)] IDA在充电时似乎对这个参数非常敏感。 请在用户脚本中使用合适的值覆盖此值（此值高度依赖于所使用的特定冷却）
% clear surface_area_per_face_for_49_layers;

%% 温度相关固相扩散的活化能 [ J / mol ]

% 正极
param.EaDps  = 5000;

% 负极
param.EaDns  = 5000;

%% 温度相关反应常数的活化能 [ J / mol ]

% 正极
param.Eakip  = 5000;

% 负极
param.Eakin  = 5000;


%% 初始条件

% 电解质锂离子初始浓度 [mol/m^3]
param.ce_init = 1000;

% 电池的初始温度 [K]
param.T_init = 298.15;

%% 模拟器参数

% 选择用于近似固相扩散的模型
% 允许值为：
% 1 - 抛物线近似（双参数模型）
%
% 2 - 高阶多项式（三参数模型）
%
% 3 - 全扩散模型
%

param.SolidPhaseDiffusion = 1;

% 选择数值方法来评估菲克扩散定律
%（即，如果 param.SolidPhaseDiffusion = 3）
% 允许值为：
% 1 - 9阶有限差分方案（至少需要10个离散点）
%
% 2 - 光谱法（未完全实现）

param.SolidPhaseDiffusionNumericalScheme = 1;

% 边界值设置/计算模式

%可以在整个代码中读取或设置 CV 边界的 % 值。 例如，输出电压的值作为正负电极边界处的 Phis 值之差获得。 
%这些值可以认为是 CV 的质心值，也可以通过插值获得。 同样的事情也适用于需要设置 Dirichlet 边界条件的情况。
%当需要设置/读取边界值时，此标志确定是插值还是考虑 CV 中心的值。

% 允许值为：
% 1 - 考虑 CV 的中心设置/计算感兴趣的值的边界值（默认值，根据相关论文）
%
% 2 - 设置/计算边界值考虑插值技术感兴趣的值
%
param.edge_values = 1;

% 集成步骤 [s]
param.sim_datalog_interval = 0.5; % 记录数据的间隔（尽管 IDA 是一个自适应求解器，在模拟结束时，用户可能会喜欢特定时间间隔的输出变量）

% 最小截止电压 [V]
param.CutoffVoltage = 2.5;

% 最大截止电压 [V]
param.CutoverVoltage = 4.3;

% 最小截止 SOC [%]
param.CutoffSOC = 0.9;

% 最大截止 SOC [%]
param.CutoverSOC = 90;

% 用于在轴向（贯穿厚度）方向离散铝集流体域的控制体积数
param.Nal   = 10; % 如果使用集总热模型则不相关

% 用于在轴向（贯穿厚度）方向离散正极域的控制体积数
param.Np    = 10;

% 用于在轴向（贯穿厚度）方向离散隔膜的控制体积数
param.Ns    = 10;

% 用于在轴向（贯穿厚度）方向离散负电极域的控制体积数
param.Nn    = 10;

% 用于在轴向（贯穿厚度）方向离散化铜集流体域的控制体积数
param.Ncu   = 10; % 如果使用集总热模型则不相关


% 如果选择全扩散模型（菲克定律），下面的两个参数定义了固体颗粒内部离散点的数量。
% 用于离散每个阴极的控制体积（壳）数
% 径向/球形方向的粒子
param.Nr_p = 10;

% 用于离散每个阳极的控制体积（壳）数
% 径向/球形方向的粒子
param.Nr_n = 10;

% 固相中锂离子的初始浓度 [mol/m^3]
param.init_cell_soc = init_cell_soc_percent/100; % 转换为 0 到 1 之间的分数
% 正极。 初始浓度自动确定
% 作为提供给 Parameters_init 脚本的 SOC 的函数
param.cs_p_init = ((param.init_cell_soc*(param.theta_max_pos-param.theta_min_pos) + param.theta_min_pos))*param.cs_maxp;

% 负极。 初始浓度自动确定
% 作为提供给 Parameters_init 脚本的 SOC 的函数
param.cs_n_init                     = ((param.init_cell_soc*(param.theta_max_neg-param.theta_min_neg) + param.theta_min_neg))*param.cs_maxn;

param.cs_neg_saturation             = (0.01*param.CutoverSOC*(param.theta_max_neg-param.theta_min_neg) + param.theta_min_neg)*param.cs_maxn;
param.enable_csneg_Saturation_limit = 0; % 此参数设置为 1 时，当负电极中任何节点的表面浓度达到饱和值时强制终止模拟（在上面的参数中设置）
param.cs_sat_thresh = 1.0;               % 快速充电算法的阈值 cs 分数

% 在 matlab 命令行中启用或禁用范围
param.Scope = 1;

% 启用或禁用标题信息的打印
param.PrintHeaderInfo = 1;

%% 外部功能

% 此字段可用作额外结构，并传递给所有外部脚本。
param.extraData = [];

% 定义必须调用以计算应用电流值的外部函数的名称。 该函数在集成过程中被调用。

param.CurrentDensityFunction    = @getInputCurrentDensity; % 包含可变电流密度曲线的句柄（外部函数文件）
param.PowerDensityFunction      = @getInputPowerDensity;     % 包含可变功率密度配置文件的句柄（外部函数文件）

% 定义用于在模拟过程中计算材料的物理和传输特性的外部函数的名称。 请参考现有函数以深入了解自定义实现。

% 电解质扩散系数
param.ElectrolyteDiffusionFunction          = @electrolyteDiffusionCoefficients;
% 电解质电导系数
param.ElectrolyteConductivityFunction       = @electrolyteConductivity;
% 开路电位
param.OpenCircuitPotentialFunction          = @openCircuitPotential;
% 固相扩散系数
param.SolidDiffusionCoefficientsFunction    = @solidPhaseDiffusionCoefficients;
% 反应率
param.ReactionRatesFunction                 = @reactionRates;

%如果提供了函数句柄，则在每个积分步骤之后调用该函数。 
%除了 extraData 结构和时序信息之外，还提供了与当前集成步骤相关的所有状态。
%以 socEstimator 为例
param.SOC_estimation_function = @socEstimator;

%% 恒电位模式
% 该值（仅在OperatingMode 标志设置为3 时适用）用于以恒电位方式控制电池。
param.V_reference = 4;

%% 公差
% 积分器 (IDA) 容差
param.AbsTol = 1e-6;
param.RelTol = 1e-6;

%% 老化参数（测试目的，测试版）

param.EnableAgeing = 0;

% 初始 SEI 电阻值 [Ohm m^2]
param.R_SEI     = 0.01;
%摩尔质量 [kg/mol]
%注意：在开发锂离子电池的第一原理容量衰减模型时，Ramadass 等人。
% M_p 的测量单位以及数字本身是错误的。 
%请参阅Review of models for predicting the cycling performance of lithium ion batterise, Santhanagopalan et al.
param.M_n               = 73e-3;
% Admittance 准入                               [S/m]
param.k_n_aging         = 3.79e-7;
% 副反应电流密度            [A/m^2]
param.i_0_jside         = 0.80e-10;
% 副反应开路电压    [V]
param.Uref_s            = 0.4;
%特定化学反应的 1C 电流 [A/m^2]
param.I1C               = 29.23;
% 老化动力学中使用的加权因子。 参见 ionic Flux.m 文件中副反应电流密度的定义。
param.w 				= 2;


% 如果用户想在计算过程中使用雅可比矩阵，则设置为 1。
param.UseJacobian       = 0;

% 该值（如果设置）表示从积分器使用的雅可比函数。 它必须是 CasADi 包的“函数”类。 
%如果提供，UseJacobian=1，它将用于加速集成过程。 
%如果未提供，当 UseJacobian=1 时，代码将自行计算雅可比行列式。
param.JacobianFunction = [];


% LIONSIMBA batteryModel.m 脚本返回的 DAE 系统类型。 此功能正在开发中。
% 允许值是：
% 1 - 方程以分析形式返回，写为隐式 DAE，即
%
% x_dot - f(x,z) | 时间微分方程 z - g(x,z) | 代数方程
%
% 2 - 方程以解析形式返回，其中时间微分方程以显式形式写入，即
%
% f(x,z) | 时间微分方程
% z - g(x,z) | 代数方程
%
param.daeFormulation = 1;


% 不要改变以下几行
% 下面几行代码用于识别LIONSIMBA执行的平台是Matlab还是Octave

% 检查代码是否在 Octave 下运行。 如果运行 Octave，下面的指令应该返回 5。如果提供 0，则表示 Matlab 正在执行
%param.isMatlab = exist('OCTAVE_VERSION', 'builtin') == 0 ;

end
