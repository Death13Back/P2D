function results = startSimulation(t0,tf,initialState,input_density,startParameters)
%	开始模拟锂离子电池.
%
%   results = startSimulation(t0,tf,initialStates,input_density,startParameters)
%
%   Input:
%       - t0 : 初始时间
%       - tf : 最终时间
%       - initialStates : 包含用于初始化的数据的结构
%       - input_density  : 施加电流/功率密度。如果是负，电池放电。如果是正，电池充电。
%       - startParameters : 如果提供的话，它定义了包含参数的单元格数组在模拟中使用的％结构。必须通过Parameters_Init脚本获得每个结构。
%		如果使用N个参数结构的单元阵列，则模拟器将串联使用N个单元进行模拟。
%		如果使用1个参数结构的单元，则将模拟单个小区。
%
%   输出:
%       - results : ％ - 结果：包含依赖的解决方案的结构变量。存储的结果与时间瞬间具有多行和列作为离散量的数量。如果多个单元格模拟，索引I用于访问第i的电池包。
%
%     results.Phis{i}:                      固相电位
%     results.Phie{i}:                      液相电位
%     results.ce{i}:                        电解质浓度
%     results.cs_surface{i}:                电极表面浓度
%     results.cs_average{i}:                电极平均浓度
%     results.time{i}:                      插值模拟时间
%     results.int_internal_time{i}:         积分器时间步长
%     results.ionic_flux{i}:                离子通量
%     results.side_reaction_flux{i}:        副反应通量
%     results.SOC{i}:                       SOC状态
%     results.SOC_estimated{i}:             根据用户定义的函数估计荷电状态
%     results.Voltage{i}:                   电池电压
%     results.Temperature{i}:               电池温度
%     results.Qrev{i}:                      可逆发热率
%     results.Qrxn{i}:                      反应生热率
%     results.Qohm{i}:                      欧姆发热率
%     results.film{i}:                      副反应膜电阻
%     results.R_int{i}:                     内阻
%     results.Up{i}:                        阴极开路电位
%     results.Un{i}:                        阳极开路电位
%     results.etap{i}:                      阴极过电位
%     results.etan{i}:                      阳极过电位
%     results.parameters{i}:                仿真模拟过程中的参数
%     results.JacobianFun:                  计算过程中的雅可比矩阵



% % % try
% % %     test_lionsimba_folder
% % % catch
% % %     error('您似乎没有将 battery_model_files 目录和其中的文件夹添加到 Matlab 路径中。 请修复此问题并重新启动模拟.')
% % % end
% % % 
% % % %
% % % 

if(isempty(startParameters))
    % 用户加载电池参数
    param{1} = Parameters_init;
else
    % 用户提供的参数
    param = startParameters;
end

% 验证输入电流密度/功率密度值
if(param{1}.OperatingMode==1 || param{1}.OperatingMode==2)
    if(~isreal(input_density) || isnan(input_density) || isinf(input_density) || isempty(input_density))
        error('The input current/power densities provided by user is complex valued or is NaN. Please check the values and restart simulation.')
    end
end

% 检查环境工具的可用性
% % % checkEnvironment(param,nargin);

% % % % 如果启用，打印标题信息
% % % if(param{1}.PrintHeaderInfo==1)
% % %     headerInfo(version)
% % % end

% 如果一切正常，开始模拟。
results = mainCore(t0,tf,initialState,input_density,param);
end

function results = mainCore(t0,tf,initialState,input_density,param)

% 存储原始参数结构，以便在模拟结束时将其返回
param_original  = param;

% 获取必须模拟的单元总数
n_cells         = length(param);

% 检查在模拟恒电位条件时是否模拟了更多电池
if(param{1}.OperatingMode==3 && n_cells~=1)
    clc;
    error('!!!ERROR!!! -- 恒电位模拟仅适用于单电池组 -- !!!ERROR!!!')
end

% 检查是否给出了初始状态结构
[Y0_existence,YP0_existence,Y0,YP0] = checkInitialStates(initialState);
%在选定的操作模式之间切换，以定义合适的 getCurrentDensity（getPowerDensity）函数。
%如果正在模拟多个电池，则从参数结构的第一个元素中检索可变电流/功率分布或其恒定值。 
%这是有效的，因为在串联串中，所有电池都取相同的电流。

switch(param{1}.OperatingMode)
    case 1
        param{1}.getCurrentDensity = @(t,t0,tf,x,param,extra)input_density;
    case 2
        param{1}.getPowerDensity = @(t,t0,tf,x,param,extra)input_density;
    case 3
        param{1}.getCurrentDensity  = @(t,t0,tf,x,param,extra)0; % 虚拟值，将被求解器覆盖
        param{1}.getPowerDensity    = @(t,t0,tf,x,param,extra)0; % 不确定这个？
    case 4
        param{1}.getCurrentDensity = param{1}.CurrentDensityFunction; % 外部，即用户提供的来自参数化文件中加载的外部文件的函数
    case 5
        param{1}.getPowerDensity = param{1}.PowerDensityFunction; % 外部，即用户提供的来自参数化文件中加载的外部文件的函数
    otherwise
        error('不支持操作模式');
end

% 定义绝对和相对容差。 如果需要更多单元格，这些值取自第一个参数结构。
opt.AbsTol      = param{1}.AbsTol;
opt.RelTol      = param{1}.RelTol;

n_diff          = zeros(n_cells,1);
n_alg           = zeros(n_cells,1);
start_x_index   = 1;
start_xp_index  = 1;

% 对于每个单元格，分配内存来存储用于估计 SOC 的外部函数。
SOC_estimate    = cell(n_cells,1);

% 对单元格进行多次检查
for i=1:n_cells
    % 检查 daeFormulation 标志以防调用 startSimulation
    if(param{i}.daeFormulation~=1)
        error(['确保将每个单元参数结构的 daeFormulation 标志设置为 1， 电池包 ', num2str(i),' 没有遵守这个限制.'])
    end

    % 当使用菲克扩散定律时，至少需要 10 个离散点。 如果不满足此条件，则引发错误。
    if((param{i}.Nr_p<10 || param{i}.Nr_n<10) && param{i}.SolidPhaseDiffusion==3)
        error('在阴极和阳极中，粒子的离散点数必须至少为 10')
    end
    % 检查是否设置了 SOC 估计功能句柄。 如果函数句柄尚未定义或没有正确数量的输入参数，则返回空值。
    if(isempty(param{i}.SOC_estimation_function) || nargin(param{i}.SOC_estimation_function)~=6)
        SOC_estimate{i} = @(a,b,c,d,e,f,g,h,i,j,k)[];
    else
        SOC_estimate{i} = @SOCestimation;
    end
    param{i}.Nsum      = param{i}.Np + param{i}.Ns + param{i}.Nn;
    param{i}.Nsum_nos  = param{i}.Np + param{i}.Nn;

    % 定义离散化步骤
    param{i}.deltax_al     = 1 / param{i}.Nal;
    param{i}.deltax_p      = 1 / param{i}.Np;
    param{i}.deltax_s      = 1 / param{i}.Ns;
    param{i}.deltax_n      = 1 / param{i}.Nn;
    param{i}.deltax_cu     = 1 / param{i}.Ncu;

    % 计算用于存储微分和代数变量位置的索引。
    param{i} = computeVariablesIndices(param{i});

    % 预分配用于固相电位的微分矩阵。 这个可以在这里完成，因为这样的矩阵被认为是时间不变的。
    param{i} = solidPhaseDifferentiationMatrices(param{i});

    % 使用电流密度/功率密度的值初始化 Param.I_density 或 Param.P_density（在上面几行的代码行中设置）
    if param{1}.OperatingMode==1 || param{1}.OperatingMode==4
        param{i}.I_density = 0;%param{1}.getCurrentDensity(0,t0,tf,x,param,param{1}.extraData);
    elseif param{1}.OperatingMode==2 || param{1}.OperatingMode==5
        param{i}.P_density = 0;%param{1}.getPowerDensity(0,t0,tf,x,param,param{1}.extraData);
    else
        param{i}.I_density = 0; % 这是一个虚拟值
        param{i}.P_density = 0;
    end

    % 预分配用于固相的微分矩阵
    % 在 Fick 定律的情况下，扩散。
    param{i} = solidPhaseDiffusionDifferentiationMatrices(param{i});

    % 获取微分状态的未知数
    [~, ~, ~, ~, ~, n_diff(i)] =   differentialInitialConditions(param{i});
    % 获取代数状态的未知数
    [~, n_alg(i), ~] = initialise_model(param{i});

    % 存储每个单元格的微分和代数变量的数量。
    param{i}.ndiff = n_diff(i);
    param{i}.nalg  = n_alg(i);

    % x_index 变量将在电池模型文件中用于索引目的。
    param{i}.x_index    = (start_x_index:n_diff(i)+n_alg(i)+start_x_index-1);

    param{i}.xp_index   = (start_xp_index:n_diff(i)+start_xp_index-1);

    % 更新（可能的）下一个单元格的起始 x_index 值
    start_x_index       = n_diff(i)+n_alg(i)+start_x_index;

    start_xp_index      = n_diff(i)+n_alg(i)+start_xp_index;
end


for i=1:n_cells
    if((Y0_existence==0) && (YP0_existence==0))
        % 获取第 i 个单元的微分状态的初始条件
        [cs_average_init, ce_init, T_init, film_init, Q_init] =   differentialInitialConditions(param{i});
        % 求解代数方程以找到代数方程的一组半一致初始条件。 这将有助于 DAE 求解器作为热启动。
        init_point = initialise_model(param{i});

        % Build the initial values array for the integrator
        Yt0 = [ce_init;cs_average_init;T_init;film_init;Q_init;init_point];
        Y0  = [Y0;Yt0];
        YP0 = [YP0;zeros(size(Yt0))];
    end
end

if(n_cells==1)
    nc = ' cell';
else
    nc = ' cells';
end

% 清空使用过的数组
ce_t            = cell(n_cells,1);
cs_bar_t        = cell(n_cells,1);
T_t             = cell(n_cells,1);
jflux_t         = cell(n_cells,1);
Phis_t          = cell(n_cells,1);
Phie_t          = cell(n_cells,1);
cs_star_t       = cell(n_cells,1);
t_tot           = cell(n_cells,1);
Qrev_t          = cell(n_cells,1);
Qrxn_t          = cell(n_cells,1);
Qohm_t          = cell(n_cells,1);
SOC_t           = cell(n_cells,1);
Voltage_t       = cell(n_cells,1);
SOC_estimated_t = cell(n_cells,1);
film_t          = cell(n_cells,1);
js_t            = cell(n_cells,1);
R_int_t         = cell(n_cells,1);
curr_density_t  = cell(n_cells,1);
Up_t            = cell(n_cells,1);
Un_t            = cell(n_cells,1);
etap_t          = cell(n_cells,1);
etan_t          = cell(n_cells,1);
dudtp_t         = cell(n_cells,1);
dudtn_t         = cell(n_cells,1);
Q_t             = cell(n_cells,1);
yp_original     = YP0';

% 该标志用于通知仿真停止。 如果为0表示一切顺利。
exit_reason     = 0;

% 定义要传递给残差函数的结构
ida_user_data.param  = param;
ida_user_data.t0     = t0;
ida_user_data.tf     = tf;

% 定义代数和微分变量。
% id:1-> 微分变量，
% id:0-> 代数变量.
id = [];
constraint_vector=[];
for i=1:n_cells
    id                                              = [id;ones(n_diff(i),1);zeros(n_alg(i),1)];
    temp_constraint_vector                          = zeros(n_diff(i)+n_alg(i),1);
    temp_constraint_vector(param{i}.Phis_indices)   = 1; % 对所有节点中的固相电势强制为正
    constraint_vector                               = [constraint_vector;temp_constraint_vector];
end

JacFun = [];

% % 该语句检查用户是否想要使用雅可比矩阵，并且（如果是）它是否已经作为参数结构的一部分提供。
% if(param{1}.UseJacobian==1 && isempty(param{1}.JacobianFunction))
%     % 如果用户想使用雅可比行列式，但参数结构中没有提供，则评估一个新的雅可比行列式。
%     disp('Evaluating the analytical form of the Jacobian matrix. Please wait...')
%     % 导入casadi框架
%     import casadi.*
%     % 定义符号变量。
%     xsym    = SX.sym('x',[sum(n_diff)+sum(n_alg),1]);
%     xpsym   = SX.sym('xp',[sum(n_diff)+sum(n_alg),1]);
%     cj      = SX.sym('cj',1);
% 
%     % 得到以符号方式以隐式形式写入的模型方程。
%     [dx_tot, ~, ~] = batteryModel(0,xsym,xpsym,ida_user_data);
% 
%     % 评估雅可比矩阵。 （请参阅sundials指南关于雅可比结构的更多信息）。
%     J = jacobian(dx_tot,xsym) + cj*jacobian(dx_tot,xpsym);
% 
%     % 为一组给定的微分和代数变量定义雅可比计算的函数。
%     JacFun = Function('fJ',{xsym,cj},{J});
% 
%     % 将函数存储到一个结构中，这样 IDA 将使用它来评估雅可比矩阵（请参阅模拟器工具中的 jacobianFunction.m）。
%     ida_user_data.fJ = JacFun;
% 
%     % Define the options for Sundials
%     options = IDASetOptions('RelTol', opt.RelTol,...
%         'AbsTol'        , opt.AbsTol,...
%         'MaxNumSteps'   , 1500,...
%         'VariableTypes' , id,...
%         'UserData'      , ida_user_data,...
%         'JacobianFn'    , @jacobianFunction,...
%         'ConstraintTypes', constraint_vector);
% 
% elseif(param{1}.UseJacobian==1 && ~isempty(param{1}.JacobianFunction))
%     % 如果需要使用雅可比，并且在参数结构中也已经提供了，直接使用即可。
% 
%     disp('用户提供的雅可比矩阵解析函数.')
%     % 如果用户提供了雅可比行列式，使用它.
%     JacFun = param{1}.JacobianFunction;
% 
%     % 将此函数指针传递给 IDA 将调用的用于评估雅可比值的例程。
%     ida_user_data.fJ = JacFun;
% 
%     % 定义Sundials选项
%     options = IDASetOptions('RelTol', opt.RelTol,...
%         'AbsTol'        , opt.AbsTol,...
%         'MaxNumSteps'   , 1500,...
%         'VariableTypes' , id,...
%         'UserData'      , ida_user_data,...
%         'JacobianFn'    , @jacobianFunction,...
%         'ConstraintTypes', constraint_vector);
% else
%     % 在这种情况下，用户不想使用雅可比矩阵。 取而代之的是计算数值近似值。
% 
%     % 定义Sundials选项
%     options = IDASetOptions('RelTol', opt.RelTol,...
%         'AbsTol'        , opt.AbsTol,...
%         'MaxNumSteps'   , 1500,...
%         'VariableTypes' , id,...
%         'UserData'      , ida_user_data,...
%         'ConstraintTypes', constraint_vector);
% end





    % 定义Sundials选项
    options = IDASetOptions('RelTol', opt.RelTol,...
        'AbsTol'        , opt.AbsTol,...
        'MaxNumSteps'   , 1500,...
        'VariableTypes' , id,...
        'UserData'      , ida_user_data,...
        'ConstraintTypes', constraint_vector);

% 仅在使用自定义当前配置文件时启用处理分段输入的功能
if(param{1}.OperatingMode == 4)
    % 添加到 IDA 选项结构指向用于识别事件存在的函数的指针（在这种特殊情况下应用的输入电流中的不连续性）
    options.RootsFn     = @rootFinder;
    options.NumRoots    = 1;
end

% 初始化求解器
IDAInit(@batteryModel,t0,Y0,YP0,options);

if param{i}.Scope == 1
    disp(['Finding a set of consistent ICs for ',num2str(n_cells),nc,' battery pack. Please wait..'])
end
% 找到一致的初始条件
[~, yy, ~] = IDACalcIC(t0+10,'FindAlgebraic');

% 初始化开始积分时间
t = t0;

% 在结果中存储初始状态值。
y = yy';

% 存储完整状态集的总数
y_tot = y;

[ce_t,cs_bar_t,T_t,jflux_t,Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t,Up_t,Un_t,R_int_t,curr_density_t,Voltage_t,SOC_estimated_t,Qrev_t,Qrxn_t,Qohm_t,Q_t,~,dudtp_t, dudtn_t,t_tot] =...
    storeSimulationResults(n_cells,ce_t,cs_bar_t,T_t,jflux_t,Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t,curr_density_t,Voltage_t,SOC_estimated_t,Up_t,Un_t,R_int_t,Qrev_t,Qrxn_t,Qohm_t,Q_t,dudtp_t, dudtn_t, t_tot, y, t,SOC_estimate,t0,tf, param);

sim_time = 0;  % 模拟时间（即wall）仅用于报告目的。 不用于控制求解器或时间循环
% 循环直到积分时间达到 tf.
while(t<tf)

    % 检查仿真停止条件是否满足
    exit_reason = checkSimulationStopConditions(n_cells, Phis_t, cs_bar_t, T_t, param);
    % 如果达到停止条件，则停止仿真
    if(exit_reason~=0)
        break;
    end

    %% 求解 DAE 集
    % 求解器 IDA 用于求解生成的 DAE 集。 有关语法及其用法的更多信息，请参阅 IDA 手册。
    tic
    [status, t, y]   = IDASolve(tf,'OneStep');

    % 如果状态 > 0，则表示在求解方程期间已找到根
    if(status>0)
        % 在自定义电流分布的情况下，根由施加的电流密度中存在的不连续性决定
        if(param{1}.OperatingMode==4)
            % 定义一个新的时刻，在该时刻使用状态的最后已知值重新初始化求解器
            t = t*(1+1e-5);
            IDAReInit(t,y,0*y,options);

            % 从新点开始寻找一致的初始条件并继续积分
            [~, y, yp] = IDACalcIC(t+10,'FindAlgebraic');
        else
            error(e.message);
        end
    end
    sim_time    = sim_time+toc;
    y           = y';
        % 存储包含所有积分的所有状态的矩阵
    % 时间步骤。
    y_tot       = [y_tot;y];
    if(status==0)
        % 在每个时间步存储导数信息
        yp_original = [yp_original;IDAGet('DerivSolution',t,1)'];
    else
        yp_original = [yp_original;yp'];
    end

    [ce_t,cs_bar_t,T_t,jflux_t,Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t,Up_t,Un_t,R_int_t,curr_density_t,Voltage_t,SOC_estimated_t,Qrev_t,Qrxn_t,Qohm_t,Q_t,tot_voltage,dudtp_t, dudtn_t,t_tot] =...
        storeSimulationResults(n_cells,ce_t,cs_bar_t,T_t,jflux_t,Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t,curr_density_t,Voltage_t,SOC_estimated_t,Up_t,Un_t,R_int_t,Qrev_t,Qrxn_t,Qohm_t,Q_t,dudtp_t, dudtn_t, t_tot, y, t,SOC_estimate,t0,tf, param);

    % 如果输出范围处于活动状态，则向用户显示其他信息
    if(param{1}.Scope==1)
        if(n_cells==1)
            temperature     = T_t{1}(end,end);
            % 如果使用菲克扩散定律，在评估 SOC 之前，需要计算每个颗粒中的平均固体浓度。
            Sout = internalSOCestimate(cs_bar_t,param,1);
            clc
            fprintf(['No. of cells in the pack \t',num2str(n_cells),'\n']);
            fprintf(['Time \t\t\t\t\t',num2str(t),' s\n']);
            % 如果恒电位模式正在运行，施加的电流将作为 DAE 的解决方案。 否则由用户提供。
            if(param{1}.OperatingMode==1) || (param{1}.OperatingMode==4)
                fprintf(['Applied current density \t\t',num2str(y(end)),' A/m^2\n']);
            elseif (param{1}.OperatingMode==2) || (param{1}.OperatingMode==5)
                fprintf(['Applied power density \t\t',num2str(param{1}.getPowerDensity(t,t0,tf,param{1}.extraData)),' W/m^2\n']);
            end

			switch param{1}.edge_values
				% 线性插值
				case 2
					output_voltage = (1.5*Phis_t{1}(end,1)-0.5*Phis_t{1}(end,2)) - (1.5*Phis_t{1}(end,end) - 0.5*Phis_t{1}(end,end-1));
				% 质心处的值
				otherwise
					output_voltage = Phis_t{1}(end,1) - Phis_t{1}(end,end);
			end
            fprintf(['Voltage \t\t\t\t',          num2str(output_voltage),   ' V\n']);
            fprintf(['Temperature \t\t\t',        num2str(temperature),                           ' K\n']);
            fprintf(['SOC \t\t\t\t\t',            num2str(Sout),                                  ' %% \n']);
            fprintf(['Cutoff Voltage \t\t\t',     num2str(param{1}.CutoffVoltage),                ' V\n']);
            fprintf(['Cutover Voltage \t\t',      num2str(param{1}.CutoverVoltage),               ' V\n']);
            fprintf(['Internal Resistance \t',    num2str(R_int_t{1}(end)),                       ' Ohm m^2\n']);
            fprintf(['Absolute tolerance \t\t',   num2str(param{1}.AbsTol),                       '\n']);
            fprintf(['Relative tolerance \t\t',   num2str(param{1}.RelTol),                       '\n']);
            fprintf(['Initial int. time \t\t',    num2str(t0),                                    ' s\n']);
            fprintf(['Final int. time \t\t',      num2str(tf),                                    ' s\n']);
            fprintf(['N. of unknowns \t\t\t',     num2str(length(y)),                             ' \n']);
        else
            clc
            fprintf(['No. of cells in the pack \t',num2str(n_cells),'\n']);
            fprintf(['Time \t\t\t\t\t',num2str(t),' s\n']);
            % 如果恒电位模式正在运行，施加的电流将作为 DAE 的解决方案。 否则由用户提供。
            if(param{1}.OperatingMode==1) || (param{1}.OperatingMode==4)
                fprintf(['Applied current density \t\t',num2str(y(end)),' A/m^2\n']);
            elseif (param{1}.OperatingMode==2) || (param{1}.OperatingMode==5)
                fprintf(['Applied power density \t\t',num2str(param{1}.getPowerDensity(t,t0,tf,param{1}.extraData)),' W/m^2\n']);
            end

            fprintf(['Voltage \t\t\t\t',          num2str(tot_voltage),       ' V\n']);
            fprintf(['Absolute tolerance \t\t',   num2str(param{1}.AbsTol),   ' \n']);
            fprintf(['Relative tolerance \t\t',   num2str(param{1}.RelTol),   ' \n']);
            fprintf(['Initial int. time \t\t',    num2str(t0),                ' s\n']);
            fprintf(['Final int. time \t\t',      num2str(tf),                ' s\n']);
            fprintf(['N. of unknowns \t\t\t',     num2str(length(y)),         ' \n']);
        end
    end
end

disp(['Elasped time: ',num2str(sim_time),' s']);

% 为固定时间步长值进行插值
t_tot_original = t_tot;

% 构建用于插值的时间向量
time_vector = (t0:param{i}.sim_datalog_interval:tf);

% 如果模拟在用户设置的最终时间之前停止，请更改 tf 变量为了插入可用值。

if(t<tf)
    % 设置最终时间等于最后一个积分步骤
    tf          = t;
    % 重新定义用于插值的时间向量
    time_vector = (t0:param{i}.sim_datalog_interval:tf);
end

% 如果至少完成了一个积分步骤，则检索一阶时间导数信息。 否则使用初始数据。
if(time_vector(end)>t0)
    % 检索最后一个时间步的导数信息
    yp          = interp1(t_tot{i},yp_original,time_vector(end))';
else
    % 如果 SUNDIALS 执行的积分步长小于参数化步长，则将初始数据作为初始状态集返回。
    yp          = YP0;
    yp_original = YP0;
end

% IDA 求解器分配的空闲内存
IDAFree

% 这些变量将用于存储积分过程的原始结果。
Phis_t_o            = cell(n_cells,1);
Phie_t_o            = cell(n_cells,1);
ce_t_o              = cell(n_cells,1);
cs_star_t_o         = cell(n_cells,1);
cs_average_t_o      = cell(n_cells,1);
jflux_t_o           = cell(n_cells,1);
SOC_t_o             = cell(n_cells,1);
T_t_o               = cell(n_cells,1);
Voltage_t_o         = cell(n_cells,1);
SOC_estimated_t_o   = cell(n_cells,1);
film_t_o            = cell(n_cells,1);
js_t_o              = cell(n_cells,1);
R_int_t_o           = cell(n_cells,1);
Up_t_o              = cell(n_cells,1);
Un_t_o              = cell(n_cells,1);
dudtp_t_o           = cell(n_cells,1);
dudtn_t_o           = cell(n_cells,1);
Qrev_t_o            = cell(n_cells,1);
Qrxn_t_o            = cell(n_cells,1);
Qohm_t_o            = cell(n_cells,1);
etap_t_o            = cell(n_cells,1);
etan_t_o            = cell(n_cells,1);
app_current_t_o     = cell(n_cells,1);
Q_t_o               = cell(n_cells,1);
y                   = [];
y_original          = [];
for i=1:n_cells
    % 保存过电位
    etap_t{i} = Phis_t{i}(:,1:param{i}.Np)-Phie_t{i}(:,1:param{i}.Np)-Up_t{i};
    if(param{i}.EnableAgeing==1)
        etan_t{i} = Phis_t{i}(:,param{i}.Np+1:end)-Phie_t{i}(:,param{i}.Np+param{i}.Ns+1:end)-Un_t{i} - param{1}.F*jflux_t{i}(:,param{i}.Np+1:end).*(param{i}.R_SEI+film_t{i}./param{i}.k_n_aging);
    else
        etan_t{i} = Phis_t{i}(:,param{i}.Np+1:end)-Phie_t{i}(:,param{i}.Np+param{i}.Ns+1:end)-Un_t{i};
    end
    if(param{i}.sim_datalog_interval>0)
        % 存储原始结果
        Phis_t_o{i}            = Phis_t{i};
        Phie_t_o{i}            = Phie_t{i};
        ce_t_o{i}              = ce_t{i};
        cs_star_t_o{i}         = cs_star_t{i};
        cs_average_t_o{i}      = cs_bar_t{i};
        jflux_t_o{i}           = jflux_t{i};
        SOC_t_o{i}             = SOC_t{i};
        T_t_o{i}               = T_t{i};
        Voltage_t_o{i}         = Voltage_t{i};
        SOC_estimated_t_o{i}   = SOC_estimated_t{i};
        film_t_o{i}            = film_t{i};
        js_t_o{i}              = js_t{i};
        R_int_t_o{i}           = R_int_t{i};
        app_current_t_o{i}     = curr_density_t{i};
        Up_t_o{i}              = Up_t{i};
        Un_t_o{i}              = Un_t{i};
        dudtp_t_o{i}           = dudtp_t{i};
        dudtn_t_o{i}           = dudtn_t{i};
        etap_t_o{i}            = etap_t{i};
        etan_t_o{i}            = etan_t{i};
        Qrev_t_o{i}            = Qrev_t{i};
        Qrxn_t_o{i}            = Qrxn_t{i};
        Qohm_t_o{i}            = Qohm_t{i};
        Q_t_o{i}               = Q_t{i};

        if(time_vector(end)>t0)
            % 对结果进行插值
            Phis_t{i}          = interp1(t_tot{i},Phis_t{i},time_vector');
            Phie_t{i}          = interp1(t_tot{i},Phie_t{i},time_vector');
            ce_t{i}            = interp1(t_tot{i},ce_t{i},time_vector');
            cs_star_t{i}       = interp1(t_tot{i},cs_star_t{i},time_vector');
            cs_bar_t{i}        = interp1(t_tot{i},cs_bar_t{i},time_vector');
            jflux_t{i}         = interp1(t_tot{i},jflux_t{i},time_vector');
            SOC_t{i}           = interp1(t_tot{i},SOC_t{i},time_vector');
            SOC_estimated_t{i} = interp1(t_tot{i},SOC_estimated_t{i},time_vector');
            Voltage_t{i}       = interp1(t_tot{i},Voltage_t{i},time_vector');
            film_t{i}          = interp1(t_tot{i},film_t{i},time_vector');
            js_t{i}            = interp1(t_tot{i},js_t{i},time_vector');
            R_int_t{i}         = interp1(t_tot{i},R_int_t{i},time_vector');
            T_t{i}             = interp1(t_tot{i},T_t{i},time_vector');
            curr_density_t{i}  = interp1(t_tot{i},curr_density_t{i},time_vector');
            Up_t{i}            = interp1(t_tot{i},Up_t{i},time_vector');
            Un_t{i}            = interp1(t_tot{i},Un_t{i},time_vector');
            Qrev_t{i}          = interp1(t_tot{i},Qrev_t{i},time_vector');
            Qrxn_t{i}          = interp1(t_tot{i},Qrxn_t{i},time_vector');
            Qohm_t{i}          = interp1(t_tot{i},Qohm_t{i},time_vector');
            etap_t{i}          = interp1(t_tot{i},etap_t{i},time_vector');
            etan_t{i}          = interp1(t_tot{i},etan_t{i},time_vector');
            dudtp_t{i}         = interp1(t_tot{i},dudtp_t{i},time_vector');
            dudtn_t{i}         = interp1(t_tot{i},dudtn_t{i},time_vector');
            Q_t{i}             = interp1(t_tot{i},Q_t{i},time_vector');
            t_tot{i}           = time_vector';
        end
    end
    % 存储结果。 如果启用了积分步骤，则存储插值数据。
    results.Phis{i}                        = Phis_t{i};
    results.Phie{i}                        = Phie_t{i};
    results.ce{i}                          = ce_t{i};
    results.cs_surface{i}                  = cs_star_t{i};
    results.cs_average{i}                  = cs_bar_t{i};
    results.time{i}                        = t_tot{i};
    results.int_internal_time{i}           = t_tot_original{i};
    results.ionic_flux{i}                  = jflux_t{i};
    results.side_reaction_flux{i}          = js_t{i};
    results.SOC{i}                         = SOC_t{i};
    results.SOC_estimated{i}               = SOC_estimated_t{i};
    results.Voltage{i}                     = Voltage_t{i};
    results.Temperature{i}                 = T_t{i};
    results.Qrev{i}                        = Qrev_t{i};
    results.Qrxn{i}                        = Qrxn_t{i};
    results.Qohm{i}                        = Qohm_t{i};
    results.film{i}                        = film_t{i};
    results.R_int{i}                       = R_int_t{i};
    results.Up{i}                          = Up_t{i};
    results.Un{i}                          = Un_t{i};
    results.etap{i}                        = etap_t{i};
    results.etan{i}                        = etan_t{i};
    results.dudtp{i}                       = dudtp_t{i};
    results.dudtn{i}                       = dudtn_t{i};
    results.Q{i}                           = Q_t{i};
    results.parameters{i}                  = param{i};
    results.JacobianFun                    = JacFun;

    % Store original data.
    results.original.Phis{i}               = Phis_t_o{i};
    results.original.Phie{i}               = Phie_t_o{i};
    results.original.ce{i}                 = ce_t_o{i};
    results.original.cs_surface{i}         = cs_star_t_o{i};
    results.original.cs_average{i}         = cs_average_t_o{i};
    results.original.ionic_flux{i}         = jflux_t_o{i};
    results.original.side_reaction_flux{i} = js_t_o{i};
    results.original.SOC{i}                = SOC_t_o{i};
    results.original.SOC_estimated{i}      = SOC_estimated_t_o{i};
    results.original.Voltage{i}            = Voltage_t_o{i};
    results.original.Temperature{i}        = T_t_o{i};
    results.original.film{i}               = film_t_o{i};
    results.original.R_int{i}              = R_int_t_o{i};
    results.original.Up{i}                 = Up_t_o{i};
    results.original.Un{i}                 = Un_t_o{i};
    results.original.etap{i}               = etap_t_o{i};
    results.original.etan{i}               = etan_t_o{i};
    results.original.Q{i}                  = Q_t_o{i};
    results.original.parameters{i}         = param_original{i};

    % 存储初始状态数据
    y           = [y;ce_t{i}(end,:)';cs_bar_t{i}(end,:)';T_t{i}(end,:)';film_t{i}(end,:)';Q_t{i}(end,:)';jflux_t{i}(end,:)';Phis_t{i}(end,:)';Phie_t{i}(end,:)';js_t{i}(end,:)';curr_density_t{i}(end)];
    % 存储初始状态原始数据
    y_original  = [y_original;ce_t_o{i}(end,:)';cs_average_t_o{i}(end,:)';T_t_o{i}(end,:)';film_t_o{i}(end,:)';Q_t_o{i}(end,:)';jflux_t_o{i}(end,:)';Phis_t_o{i}(end,:)';Phie_t_o{i}(end,:)';js_t_o{i}(end,:)';app_current_t_o{i}(end)];
end

% 存储最后结果的数组
results.Y                           = y;
results.YP                          = yp;

results.original.Y                  = y_tot;
results.original.YP                 = yp_original;

results.original.initialState.Y     = y_original;
results.original.initialState.YP    = yp_original(end,:);

results.initialState.Y              = y;
results.initialState.YP             = yp;

% 存储模拟时间
results.simulation_time             = sim_time;

%退出原因
results.exit_reason                 = exit_reason;

% 检查操作模式并相应地存储结果。
if(param{1}.OperatingMode==2)
    results.power_density          = param{1}.P_density * ones(size(t_tot{1},1),1);
    results.original.power_density = param{1}.P_density * ones(size(t_tot_original{1},1),1);
    % 可变电流曲线
elseif(param{1}.OperatingMode==5)
    results.power_density          = param{1}.getPowerDensity(t_tot{1},t0,tf,param{1}.extraData);
    results.original.power_density = param{1}.getPowerDensity(t_tot_original{1},t0,tf,param{1}.extraData);
end

results.curr_density          = curr_density_t{1};
results.original.curr_density = app_current_t_o{1};
end

function estimate = SOCestimation(t,t0,tf,param,ce_t,cs_bar_t,cs_star_t,Phie_t,Phis_t,jflux_t,T_t)
% 构建传递给函数的状态结构
states.ce           = ce_t(end,:);
states.cs_average   = cs_bar_t(end,:);
states.cs_surface   = cs_star_t(end,:);
states.Phie         = Phie_t(end,:);
states.Phis         = Phis_t(end,:);
states.ionic_flux   = jflux_t(end,:);
states.Temperature  = T_t(end,:);
% 调用估计过程
estimate = param.SOC_estimation_function(t,t0,tf,states,param.extraData,param);
end
