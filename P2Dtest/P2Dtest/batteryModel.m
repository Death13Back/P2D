function [dx_tot, flag, new_data] = batteryModel(t,x_tot,xp_tot,ida_user_data)
% BatteryModel 返回完整电池模型的残差集。



% 这两个标志未使用，但 IDA(s) 求解器需要。
flag            = 0;
new_data        = [];
% 从 IDA(s) 的 UserData 字段中检索数据
param_tot       = ida_user_data.param;
t0              = ida_user_data.t0;
tf              = ida_user_data.tf;
% 清空残差的总数组。
dx_tot          = [];

% 查找堆栈中的单元格数量
n_cells = length(param_tot);

% 评估所有电池。
for i=1:n_cells
    % 将与当前单元格相关的参数结构关联到“param”变量。
    param       = param_tot{i};
    % 从变量总数中，提取与当前单元格相关的变量。
    x       = x_tot(param.x_index);
    xp      = xp_tot(param.xp_index);

    % 所有串联的电池都承载相同的电流。
     % 尽管如此，每个单元都可以用自己的电流驱动。 这对于需要电荷平衡的 ABMS 开发很有用。

     % 检索微分变量
     % 电解质浓度
    ce          = x(param.ce_indices);
    % 平均固相浓度
    cs_barrato  = x(param.cs_average_indices);
    % 温度
    T           = x(param.T_indices);
    % 薄膜电阻
    film        = x(param.film_indices);
    % 平均通量
    Q           = x(param.Q_indices);

    % 离子通量
    jflux       = x(param.jflux_indices);
    % 固相电位
    Phis        = x(param.Phis_indices);
    % 电解质电位
    Phie        = x(param.Phie_indices);
    % 副反应通量
    js          = x(param.js_indices);
    % 当前密度
    I_density   = x(param.curr_dens_indices);

    % 电解质浓度导数
    dCe         = xp(param.ce_indices,1);
    % 平均表面浓度导数
    dCs         = xp(param.cs_average_indices,1);
    % 温度导数
    dT          = xp(param.T_indices,1);
    % 膜厚
    dFilm       = xp(param.film_indices,1);
    % 平均通量
    dQ          = xp(param.Q_indices,1);

    if param.OperatingMode==1 || param.OperatingMode==4
        param.I_density = param_tot{1}.getCurrentDensity(t,t0,tf,x_tot,param_tot,param_tot{i}.extraData); %param.I_density 通过这一重要代码行的每个时间步长都会更新，即在时间 't' 获取值（常量或来自外部函数文件）
    elseif param.OperatingMode==2 || param.OperatingMode==5
        param.P_density = param_tot{1}.getPowerDensity(t,t0,tf,param_tot{i}.extraData);   %param.P_density 通过这一重要代码行的每个时间步长都会更新，即在时间 't' 获取值（常量或来自外部函数文件）
    elseif param.OperatingMode==3    % 检查是否在恒电位充电模式下运行
        param.I_density = I_density;
    else
        error('Not a valid operating mode');
    end

    % 检查电池数量，以便为热动力学提供正确的 BC（如果热动力学处于活动状态）。 
    %当热动力学开启时，两个相邻单元之间的热通量必须保持连续。
%     if(param.TemperatureEnabled >= 1)
%         if(n_cells==1)
%             % 如果它只存在一个电池，那么该电池将与周围环境进行热交换。
%             param.T_BC_sx = param.hcell*(param.Tref-T(1))/(param.deltax_al*param.len_al);
%             param.T_BC_dx = -param.hcell*(T(end)-param.Tref)/(param.deltax_cu*param.len_cu);
%         else
%             % 如果使用多个电池，则热交换将相应地与相邻电池进行。 
%             %只有外部阴极和阳极会随着周围环境散热，而内部电池将保持热流连续。
%             if(i==1)
%                 param.T_BC_sx = param.hcell*(param.Tref-T(1))/(param.deltax_al*param.len_al);
% 
%                 % 找到当前电池最后一个体积的deltax
%                  deltax_sx = param.len_cu*param.deltax_cu;
%                  % 找到下一个电池的第一个体积的 deltax
%                  deltax_dx = param_tot{i+1}.len_al*param_tot{i+1}.deltax_al;
%                  % 求两个不同电池中两个相邻节点之间的总距离。
%                 deltax_tot = (deltax_sx/2+deltax_dx/2);
%                 % 求传热系数的调和平均值
%                 beta = (deltax_sx/2)/deltax_tot;
%                 Lambda_star = (param.Lambda_cu*param_tot{i+1}.Lambda_al)/(beta*param_tot{i+1}.Lambda_al + (1-beta)*param.Lambda_cu);
%                 % 在当前单元格的右侧设置 BC。 注意 x_next_cell(param_tot{i+1}.T_indices(1)) 指的是下一个电池的 T(1)
%                 x_next_cell   = x_tot(param_tot{i+1}.x_index);
%                 param.T_BC_dx = Lambda_star*((x_next_cell(param_tot{i+1}.T_indices(1))-T(end))/deltax_tot)/(param.deltax_cu*param.len_cu);
%             elseif(i>1 && i<n_cells)
%                 % 查找前一个电池的最后一个体积的 deltax
%                 deltax_sx = param_tot{i-1}.len_cu*param_tot{i-1}.deltax_cu;
%                 % 查找当前电池的第一个体积的 deltax
%                 deltax_dx = param.len_al*param.deltax_al;
%                 % 求两个不同单元格中两个相邻节点之间的总距离。
%                 deltax_tot = (deltax_sx/2+deltax_dx/2);
%                 % 求传热系数的调和平均值
%                 beta = (deltax_sx/2)/deltax_tot;
%                 Lambda_star = (param.Lambda_al*param_tot{i-1}.Lambda_cu)/(beta*param.Lambda_al + (1-beta)*param_tot{i-1}.Lambda_cu);
%                 % 设置当前电池左侧的BC
%                 x_prev_cell   = x_tot(param_tot{i-1}.x_index);
%                 param.T_BC_sx = -Lambda_star*((T(1)-x_prev_cell(param_tot{i-1}.T_indices(end)))/deltax_tot)/(param.deltax_al*param.len_al);
% 
%                 % 查找当前电池的最后一个体积的 deltax
%                 deltax_sx = param.len_cu*param.deltax_cu;
%                 % 找到下一个电池的第一个体积的 deltax
%                 deltax_dx = param_tot{i+1}.len_al*param_tot{i+1}.deltax_al;
%                 % 求两个不同单元格中两个相邻节点之间的总距离。
%                 deltax_tot = (deltax_sx/2+deltax_dx/2);
%                 % 求传热系数的调和平均值
%                 beta = (deltax_sx/2)/deltax_tot;
%                 Lambda_star = (param.Lambda_cu*param_tot{i+1}.Lambda_al)/(beta*param_tot{i+1}.Lambda_al + (1-beta)*param.Lambda_cu);
%                 % 设置当前电池右侧的BC
%                 x_next_cell   = x_tot(param_tot{i+1}.x_index);
%                 param.T_BC_dx = Lambda_star*((x_next_cell(param_tot{i+1}.T_indices(1))-T(end))/deltax_tot)/(param.deltax_cu*param.len_cu);
%             else
%                 % 查找前一个电池的最后一个体积的 deltax
%                 deltax_sx = param_tot{i-1}.len_cu*param_tot{i-1}.deltax_cu;
%                 % 查找当前电池的第一个体积的 deltax
%                 deltax_dx = param.len_al*param.deltax_al;
%                 % 求两个不同单元格中两个相邻节点之间的总距离。
%                 deltax_tot = (deltax_sx/2+deltax_dx/2);
%                 % 求传热系数的调和平均值
%                 beta = (deltax_sx/2)/deltax_tot;
%                 Lambda_star = (param.Lambda_al*param_tot{i-1}.Lambda_cu)/(beta*param.Lambda_al + (1-beta)*param_tot{i-1}.Lambda_cu);
%                 % 设置当前电池左侧的BC
%                 x_prev_cell   = x_tot(param_tot{i-1}.x_index);
%                 param.T_BC_sx = -Lambda_star*((T(1)-x_prev_cell(param_tot{i-1}.T_indices(end)))/deltax_tot)/(param.deltax_al*param.len_al);
%                 param.T_BC_dx = -param.hcell*(T(end)-param.Tref)/(param.deltax_cu*param.len_cu);
%             end
%         end
%     end

    % 为代数方程构建初始条件数组
    x0_alg      = [jflux;Phis;Phie;js;I_density];

    % 从代数方程中获取残差
    [dxalg,Up,Un,dudt_p,dudt_n,Keff,J_S] = algebraicStates(x0_alg,ce,cs_barrato,Q,T,film,param);

    % 评估电解质浓度的残差
    [resCe,rhsCe] = electrolyteDiffusion(ce,dCe,jflux + [zeros(param.Np,1);J_S],T,param);
    % 评估平均表面浓度的残差
    [resCs, rhsCs] = electrodeConcentration(dCs,cs_barrato,T,jflux,param);

    % 检查要使用的固体扩散近似的类型。
    if(param.SolidPhaseDiffusion==2)
        [resdQ, rhsQ] = volumeAveragedConcentrationFlux(dQ,Q,jflux,T,param);
    else
        resdQ   = dQ;
        rhsQ    = zeros(length(dQ),1);
    end

    % 老化效应仅在充电过程中发生（如此处假设）。 需要在施加的电流密度是数字量的情况和符号量的情况之间切换。
    if(isa(I_density,'casadi.SX') && param.EnableAgeing==1)
        [resDfilm, rhsDfilm] = SEI_layer(dFilm, J_S, param);         % 膜厚

        % 使用 if_else CasADi 语句根据施加的电流密度的符号确定薄膜厚度动态的切换行为。
        % resDfilm = if_else(I_density>=0,resDfilm,zeros(size(resDfilm,1),1));
        rhsDfilm = if_else(I_density>=0,rhsDfilm,zeros(size(rhsDfilm,1),1));
    elseif(param.EnableAgeing==1 && I_density >0) % 如果施加的电流是数字量，则执行常规计算

        [resDfilm, rhsDfilm] = SEI_layer(dFilm, J_S, param);   % 膜厚
    else
        resDfilm    = dFilm;
        rhsDfilm    = zeros(length(resDfilm),1);
    end
        resdT   = dT;
        rhsT    = zeros(length(dT),1);
    % 检查是否启用了热动力学。
%     switch param.TemperatureEnabled
% 	case 1 % 基于 PDE 的热模型
%          % 评估温度的残差
%         [resdT, rhsT] = thermalModel_pde(ce,Phie,Phis,Keff,jflux+ [zeros(param.Np,1);J_S],T,dT,Up,Un,dudt_p,dudt_n,x_tot(end),param);
% 	case 2 % 降阶集总热模型
%            % 表面浓度
%         cs_star = surfaceConcentration(cs_barrato,jflux,Q,T,param); % 列向量
%         cs_star_avg_pos = cs_star(1:param.Np)'*ones(param.Np,1)/param.Np;        % （标量）正极中的平均表面浓度
%         cs_star_avg_neg = cs_star(param.Np+1:end)'*ones(param.Nn,1)/param.Nn;    % （标量）负极中的平均表面浓度
%         cs_star_avg = [cs_star_avg_pos*ones(param.Np,1);cs_star_avg_neg*ones(param.Nn,1)]; % 适当复制平均表面浓度并连接
% 
%         % 计算温度的残差向量
%         [resdT, rhsT] = thermalModel_lumped(ce,cs_star_avg,Phis,Keff,jflux + [zeros(param.Np,1);J_S],T,dT,Up,Un,dudt_p,dudt_n,param,I_density);
% 	otherwise % 热动力学禁用
%         % 返回温度的常数值
%         resdT   = dT;
%         rhsT    = zeros(length(dT),1);
%     end

    if(param_tot{1}.daeFormulation==1)
        % 返回残差数组
        dx_tot = [dx_tot;resCe;resCs;resdT;resDfilm;resdQ;dxalg];
    elseif(param_tot{1}.daeFormulation==2)
        % 返回方程的 RHS
        dx_tot = [dx_tot;rhsCe;rhsCs;rhsT;rhsDfilm;rhsQ;dxalg];
    end
end
end
