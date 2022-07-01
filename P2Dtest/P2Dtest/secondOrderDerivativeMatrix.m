function [derivativeMatrix,r12dxs,dx]=secondOrderDerivativeMatrix(xl,xu,n)
% secondOrderDerivativeMatrix 使用 6 个点预先计算用于实现数值微分的矩阵。



%  网格间距
dx 		= (xu-xl)/(n-1);
%
r12dxs 	= 1./(12.0*dx^2);

% 定义用于构建数值微分矩阵的块。
mid_block        = zeros(n-4,n);

%% 数值算法
first_row       = [-415/6   +96     -36     +32/3       -3/2    0];
second_row      = [+10      -15     -4      +14         -6      +1];
i_th_row        = [-1       +16     -30     +16         -1];
semi_last_row   = [+1       -6      +14     -4          -15     +10];
last_row        = [0    -3/2    +32/3   -36     +96     -415/6];

%% 块构建
first_block = [first_row;second_row];
first_block = [first_block zeros(2,n-6)];


last_block 	= [semi_last_row;last_row];
last_block 	= [zeros(2,n-6) last_block];

row_index 	= 1;
for i=3:n-2
    mid_block(row_index,row_index:row_index+4) 	= i_th_row;
    row_index 									= row_index+1;
end

% 构建整体矩阵。 乘以 6 是为了避免由于原始矩阵定义中存在的有理值引起的数值问题。
% 在 FDM9orderElectrodeDiffusion.m 中除以 6 以获得适当的二阶导数
derivativeMatrix = 6*[first_block;mid_block;last_block];
end