function sign_input_density = evaluate_sign_input_density(param)
% 根据操作模式，evaluate_sign_input_density 返回输入电流/功率密度的符号



if param.OperatingMode==1 || param.OperatingMode==4
    sign_input_density = sign(param.I_density);
elseif param.OperatingMode==2 || param.OperatingMode==5
    sign_input_density = sign(param.P_density);
elseif param.OperatingMode==3
    sign_input_density = 1; % CV 模式的虚拟值（因为通常我们有 CV 充电）
else
    error('Not a valid operating mode.');
end

end
