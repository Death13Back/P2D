function param = solidPhaseDiffusionDifferentiationMatrices(param)
%	solidPhaseDiffusionDifferentiationMatrices 预先计算求解固相扩散方程的数值方案。 当考虑 Fick 扩散定律时将使用这些数据。



param.Rad_position_p  = linspace(0,param.Rp_p,param.Nr_p);
if(size(param.Rad_position_p,1)<size(param.Rad_position_p,2))
    param.Rad_position_p = param.Rad_position_p';
end

param.Rad_position_n  = linspace(0,param.Rp_n,param.Nr_n);
if(size(param.Rad_position_n,1)<size(param.Rad_position_n,2))
    param.Rad_position_n = param.Rad_position_n';
end

% 预先计算用于计算一阶数值导数的矩阵。 由于数值原因，粒子的半径域 [0,R] 在 [0,1] 中归一化。 相应地实施数值公式。
[param.FO_D_p,param.FO_D_c_p] = firstOrderDerivativeMatrix(0,1,param.Nr_p);
[param.FO_D_n,param.FO_D_c_n] = firstOrderDerivativeMatrix(0,1,param.Nr_n);

% 预先计算用于计算二阶数值导数的矩阵。 
%由于数值原因，粒子的半径域 [0,R] 在 [0,1] 中归一化。 
%相应地实施数值公式。
[param.SO_D_p,param.SO_D_c_p,param.SO_D_dx_p] = secondOrderDerivativeMatrix(0,1,param.Nr_p);
[param.SO_D_n,param.SO_D_c_n,param.SO_D_dx_n] = secondOrderDerivativeMatrix(0,1,param.Nr_n);

end