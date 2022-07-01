function Sout = socEstimator(t,t0,tf,states,extraData,param)
%	socEstimator 估计电池的 SOC。
% 可以在每个积分步骤之后（在时间步进循环内）调用此函数。
% 可以通过适当更改 Parameters_init 中的函数句柄来实现自定义 SOC 估计。 
%必须遵守函数的符号并返回标量值。 Sout 值将被连接起来，并在模拟结束时在结果数组中可用。



if(param.SolidPhaseDiffusion~=3)
    cs_average = states.cs_average(end,param.Np+1:end);
else
    start_index = param.Nr_p*param.Np+1;
    end_index   = start_index+param.Nr_n-1;
    cs_average  = zeros(param.Nn,1);
    for n=1:param.Nn
        cs_average(n)   = 1/param.Rp_n*(param.Rp_n/param.Nr_n)*sum(states.cs_average(end,start_index:end_index));
        start_index     = end_index + 1;
        end_index       = end_index + param.Nr_n;
    end
end
% 根据状态的当前信息估计SOC。
Csout= sum(cs_average);
Sout = 100*(1/param.len_n*(param.len_n/(param.Nn))*Csout/param.cs_maxn);

end
