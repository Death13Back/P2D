function [dPhis, dPhie, dCe] = ThermalDerivatives(Phis,Phie,ce,param)
% ThermalDerivatives 评估热力学中使用的空间导数集。
% 此函数评估热（发热）计算的相关变量的空间导数（梯度）。 请注意，空间导数是在控制体积内而不是在控制体积的边缘处计算的。



%对于下面计算的每个数值导数，
%第一个和最后一个控制体积是用一阶精度（分别是前向和后向差分方案）评估的，
%而中间控制体积近似值使用二阶精度（中心差分方案）。

%% 固相电势导数

% 正极
dPhisp = [(-3*Phis(1)+4*Phis(2)-Phis(3))/(2*param.deltax_p*param.len_p);... 						% 向前    差分
    (Phis(3:param.Np)-Phis(1:param.Np-2)) / (2*param.deltax_p*param.len_p);...						% 中心    差分
    (3*Phis(param.Np)-4*Phis(param.Np-1)+Phis(param.Np-2)) / (2*param.deltax_p*param.len_p);...		% 向后    差分
    ];

% 负极
dPhisn = [(-3*Phis(param.Np+1)+4*Phis(param.Np+2)-Phis(param.Np+3))/(2*param.deltax_n*param.len_n);... 	% 向前    差分
    (Phis(param.Np+3:end)-Phis(param.Np+1:end-2)) / (2*param.deltax_n*param.len_n);... 					% 中心    差分
    (3*Phis(end)-4*Phis(end-1)+Phis(end-2)) / (2*param.deltax_n*param.len_n);...						% 向后    差分
    ];

dPhis = [dPhisp;dPhisn];

%% 电解质电位导数

% 正极

dPhiep = [ (-3*Phie(1)+4*Phie(2)-Phie(3))/(2*param.deltax_p*param.len_p);...		% 向前 差分
    (Phie(3:param.Np)-Phie(1:param.Np-2))/(2*param.deltax_p*param.len_p);...	  	% 中心 差分
    ];

%注意！ 正电极的最后一个体积将涉及一个体积的隔膜，用于计算导数。 因此，必须对 deltax 量进行适当的考虑。

% 正极中的最后一个 CV：具有中心差分算法的导数近似
dPhie_last_p = 2*(Phie(param.Np+1)-Phie(param.Np-1))/(3 * param.deltax_p*param.len_p + param.deltax_s*param.len_s);

% 隔膜

% 注意！ 分离器的第一个体积将涉及一个体积的正截面用于计算导数。 因此，必须对 deltax 量进行适当的考虑。

% 分离器中的第一个 CV：使用中心差分算法的导数近似
dPhie_first_s = 2*(Phie(param.Np+2)-Phie(param.Np))/(param.deltax_p*param.len_p + 3* param.deltax_s*param.len_s);

% 中心差分是否
dPhies =  (Phie(param.Np+3:param.Np+param.Ns)-Phie(param.Np+1:param.Np+param.Ns-2))/(2*param.deltax_s*param.len_s);

% 注意力！ 分离器的最后一个体积将涉及负部分的一个体积，用于计算导数。 因此，必须对 deltax 量进行适当的考虑。

% 分离器中的最后一个 CV：具有中心方案的导数近似
dPhie_last_s = 2*(Phie(param.Np+param.Ns+1)-Phie(param.Np+param.Ns-1))/( param.deltax_n*param.len_n + 3*param.deltax_s*param.len_s);

% 负极

% 注意力！ 负极的第一体积将涉及一个体积的隔板部分，用于计算导数。 因此，必须对 deltax 量进行适当的考虑。

% 负极中的第一个 CV：具有中心差分的导数近似
dPhie_first_n = 2*(Phie(param.Np+param.Ns+2)-Phie(param.Np+param.Ns))/(3 * param.deltax_n*param.len_n + param.deltax_s*param.len_s);

% 中心差分算法
dPhien = [(Phie(param.Np+param.Ns+3:end)-Phie(param.Np+param.Ns+1:end-2))/(2*param.deltax_n*param.len_n);...
    (3*Phie(end)-4*Phie(end-1)+Phie(end-2))/(2*param.deltax_n*param.len_n)
    ];

dPhie = [dPhiep;dPhie_last_p;dPhie_first_s;dPhies;dPhie_last_s;dPhie_first_n;dPhien];

%% 电解质浓度导数

% 正极

dCep = [ (-3*ce(1)+4*ce(2)-ce(3))/(2*param.deltax_p*param.len_p);... 		% 向前差分算法
    (ce(3:param.Np)-ce(1:param.Np-2))/(2*param.deltax_p*param.len_p);... 	% 中心差分算法
    ];

% 注意！ 正电极的最后一个体积将涉及一个体积的隔板，用于计算导数。 因此，必须对 deltax 量进行适当的考虑。

% 正极中的最后一个 CV：具有中心差分算法的导数近似
dCe_last_p = 2*(ce(param.Np+1)-ce(param.Np-1))/(3 * param.deltax_p*param.len_p + param.deltax_s*param.len_s);

% 隔膜

% 注意！ 分离器的第一个体积将涉及一个体积的正截面用于计算导数。 因此，必须对 deltax 量进行适当的考虑。

% 分离器中的第一个 CV：具有中心差分算法的导数近似
dCe_first_s = 2*(ce(param.Np+2)-ce(param.Np))/( param.deltax_p*param.len_p + 3* param.deltax_s*param.len_s);

% 中心差分算法
dCes = (ce(param.Np+3:param.Np+param.Ns)-ce(param.Np+1:param.Np+param.Ns-2))/(2*param.deltax_s*param.len_s);

% 注意力！ 分离器的最后一个体积将涉及负部分的一个体积，用于计算导数。 因此，必须对 deltax 量进行适当的考虑。

% 分离器中的最后一个 CV：具有中心差分的导数近似
dCe_last_s = 2*(ce(param.Np+param.Ns+1)-ce(param.Np+param.Ns-1))/( param.deltax_n*param.len_n + 3*param.deltax_s*param.len_s);

% 负极

% 注意！ 负极的第一体积将涉及一个体积的隔板部分，用于计算导数。 因此，必须对 deltax 量进行适当的考虑。

% 负极中的第一个 CV：具有中心差分算法的导数近似
dCe_first_n = 2*(ce(param.Np+param.Ns+2)-ce(param.Np+param.Ns))/(3 * param.deltax_n*param.len_n + param.deltax_s*param.len_s);

dCen = [(ce(param.Np+param.Ns+3:end)-ce(param.Np+param.Ns+1:end-2))/(2*param.deltax_p*param.len_p);... 	% 中心差分算法
    (3*ce(end)-4*ce(end-1)+ce(end-2))/(2*param.deltax_n*param.len_n) 							% 向后差分算法
    ];

dCe = [dCep;dCe_last_p;dCe_first_s;dCes;dCe_last_s;dCe_first_n;dCen];

end
