function param = solidPhaseDifferentiationMatrices(param)
%solidPhaseDifferentiationMatrices 预分配用于计算 Phis 指数的数值导数的矩阵。
%鉴于在代码中，固相电导率在整个单元长度上被认为是恒定的，离散矩阵只要组装一次

% 正极矩阵
c = ones(param.Np-1,1);
d = -2*ones(param.Np,1);

A_p 				= gallery('tridiag',c,d,c);
A_p(1,1) 			= -1;
A_p(end,end-1:end) 	= [1 -1];

%负极矩阵
c = ones(param.Nn-1,1);
d = -2*ones(param.Nn,1);

A_n 			= gallery('tridiag',c,d,c);
A_n(1,1) 		= -1;
A_n(end,end) 	= -1;
% 将矩阵存储在 param 结构中以备将来使用
%包中的单元格。
if param.OperatingMode==2 || param.OperatingMode==5 % 如果在功率密度模式下运行
    A_n(end,:)=[]; % 在 algebraicStates.m 文件中应用非线性边界条件（对于复杂的功率密度 BC 推导）
end

param.A_p = A_p;
param.A_n = A_n;
end
