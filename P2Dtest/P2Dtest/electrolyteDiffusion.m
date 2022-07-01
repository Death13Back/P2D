function [resCe, rhsCe] = electrolyteDiffusion(ce,dCe,jflux,T,param)
% 电解质扩散评估电解质溶液中锂离子电解质浓度的残留量。



% 扩散系数 出于基准测试目的对此进行评论
Deff_p = param.ElectrolyteDiffusionFunction(ce(1:param.Np),T(param.Nal+1:param.Nal+param.Np),param,'p');
Deff_s = param.ElectrolyteDiffusionFunction(ce(param.Np+1:param.Np+param.Ns),T(param.Nal+param.Np+1:param.Nal+param.Np+param.Ns),param,'s');
Deff_n = param.ElectrolyteDiffusionFunction(ce(param.Np+param.Ns+1:end),T(param.Nal+param.Np+param.Ns+1:end-param.Ncu),param,'n');

% 出于基准测试目的取消注释
% Deff_p = repmat(param.Dp*param.eps_p^param.brugg_p,param.Np,1);
% Deff_s = repmat(param.Ds*param.eps_s^param.brugg_s,param.Ns,1);
% Deff_n = repmat(param.Dn*param.eps_n^param.brugg_n,param.Nn,1);

% 扩散系数的插值
[Deff_p, Deff_s, Deff_n] = interpolateDiffusionCoefficients(Deff_p,Deff_s,Deff_n,param);

%% 正极矩阵
A_p = - diag(Deff_p);

A_p(2:end,2:end)    = A_p(2:end,2:end) - diag(Deff_p(1:end-1));
A_p(1:end-1,2:end)  = A_p(1:end-1,2:end) + diag(Deff_p(1:end-1));
A_p(2:end,1:end-1)  = A_p(2:end,1:end-1) + diag(Deff_p(1:end-1));
%% 隔膜 A 矩阵
A_s = - diag(Deff_s);

A_s(2:end,2:end)    = A_s(2:end,2:end) - diag(Deff_s(1:end-1));
A_s(1:end-1,2:end)  = A_s(1:end-1,2:end) + diag(Deff_s(1:end-1));
A_s(2:end,1:end-1)  = A_s(2:end,1:end-1) + diag(Deff_s(1:end-1));
%% 负极矩阵
A_n = - diag(Deff_n);

A_n(2:end,2:end)    = A_n(2:end,2:end) - diag(Deff_n(1:end-1));
A_n(1:end-1,2:end)  = A_n(1:end-1,2:end) + diag(Deff_n(1:end-1));
A_n(2:end,1:end-1)  = A_n(2:end,1:end-1) + diag(Deff_n(1:end-1));
% 修复 A_n 的最后一个元素
A_n(end,end-1:end)  = [Deff_n(end-1) -Deff_n(end-1)];
%% A_tot 矩阵
A_tot = blockDiagonalMatrix(param,A_p,A_s,A_n);

% 除以 deltax 和正极的长度
A_tot(1:param.Np,1:param.Np)    = A_tot(1:param.Np,1:param.Np)/(param.deltax_p^2*param.len_p^2);
% 为接口条件重置线路上的值
A_tot(param.Np,:)           = 0;
A_tot(param.Np+1,:)         = 0;

% 除以 deltax 和分隔符的长度
A_tot(param.Np+1:param.Np+param.Ns,param.Np+1:param.Np+param.Ns)    = A_tot(param.Np+1:param.Np+param.Ns,param.Np+1:param.Np+param.Ns)/(param.deltax_s^2 * param.len_s^2);
% 为接口条件重置线路上的值
A_tot(param.Np+param.Ns,:)      = 0;
A_tot(param.Np+param.Ns+1,:)    = 0;

% 除以 deltax 和负极的长度
A_tot(param.Np+param.Ns+1:end,param.Np+param.Ns+1:end)      = A_tot(param.Np+param.Ns+1:end,param.Np+param.Ns+1:end)/(param.deltax_n^2 * param.len_n^2);

%% 隔膜和正极之间的界面（正极中的最后一个体积单元）

% 计算接口处的分母
den_s   = (param.deltax_p*param.len_p/2 + param.deltax_s*param.len_s/2);
% 正极末次扩散系数
last_p  = Deff_p(end-1)/(param.deltax_p*param.len_p);
% 界面扩散系数
first_s = Deff_p(end)/den_s;
% 修复边界处的值
A_tot(param.Np,param.Np-1:param.Np+1) = [last_p -(last_p+ first_s) first_s]/(param.deltax_p*param.len_p*param.eps_p);

%% 隔膜和正极之间的界面（隔膜中的第一体积）

% 计算接口处的分母
den_s       = (param.deltax_p*param.len_p/2 + param.deltax_s*param.len_s/2);
% 隔板中的第一扩散系数
second_s    = Deff_s(1)/(param.deltax_s*param.len_s);
% 界面扩散系数
first_s     = Deff_p(end)/den_s;

A_tot(param.Np+1,param.Np:param.Np+2) = [first_s -(first_s+second_s) second_s]/(param.deltax_s*param.len_s*param.eps_s);

%% 隔膜和负极之间的界面（隔膜中的最后体积）

% 计算接口处的分母
den_s   = (param.deltax_s*param.len_s/2 + param.deltax_n*param.len_n/2);
% 隔膜中中的最后扩散系数
last_s  = Deff_s(end-1)/(param.deltax_s*param.len_s);
% 界面扩散系数
first_n = Deff_s(end)/den_s;

A_tot(param.Np+param.Ns,param.Np+param.Ns-1:param.Np+param.Ns+1) = [last_s -(last_s+first_n) first_n]/(param.deltax_s*param.len_s*param.eps_s);

%% 隔膜和负极之间的界面（负极中的第一体积单元）

% 计算接口处的分母
den_n       = (param.deltax_s*param.len_s/2 + param.deltax_n*param.len_n/2);
% 负极中的第一扩散系数
second_n    = Deff_n(1)/(param.deltax_n*param.len_n);
% 界面扩散系数
first_n     = Deff_s(end)/den_n;

A_tot(param.Np+param.Ns+1,param.Np+param.Ns:param.Np+param.Ns+2) = [first_n -(first_n+second_n) second_n]/(param.deltax_n*param.len_n*param.eps_n);

%% Useful stuff
a_tot       = [
    repmat(param.a_i(1),param.Np,1);...
    zeros(param.Ns,1);...
    repmat(param.a_i(3),param.Nn,1)...
    ];

jflux_tot   = [
    jflux(1:param.Np);...
    zeros(param.Ns,1);...
    jflux(param.Np+1:end)...
    ];

eps_tot     = [
    repmat(param.eps_p,param.Np,1);...
    repmat(param.eps_s,param.Ns,1);...
    repmat(param.eps_n,param.Nn,1)
    ];

K = 1./eps_tot;
A_eps = diag(K);

% 建立孔隙率矩阵
if(~isa(eps_tot,'casadi.MX') && ~isa(eps_tot,'casadi.SX'))
    A_eps = A_eps + diag(K(1:end-1),1);
    A_eps = A_eps + diag(K(1:end-1),-1);
    A_eps = sparse(A_eps);
else
    A_u_l                   = diag(K(1:end-1));
    A_eps(1:end-1,2:end)    = A_eps(1:end-1,2:end)+A_u_l;
    A_eps(2:end,1:end-1)    = A_eps(2:end,1:end-1)+A_u_l;
end

A_eps(param.Np,param.Np-1:param.Np+1) = 1;
A_eps(param.Np+1,param.Np:param.Np+2) = 1;

A_eps(param.Np+param.Ns,param.Np+param.Ns-1:param.Np+param.Ns+1) = 1;
A_eps(param.Np+param.Ns+1,param.Np+param.Ns:param.Np+param.Ns+2) = 1;

G = A_eps.*A_tot;

% 写出方程的 RHS
rhsCe = (G*ce + K.*(1-param.tplus).*a_tot.*jflux_tot);

% 写出方程的残差
resCe = dCe - rhsCe;
end
