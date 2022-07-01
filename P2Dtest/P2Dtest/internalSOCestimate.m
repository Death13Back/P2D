function Sout = internalSOCestimate(cs_average_t,param,i)
%   internalSOCestimate 用于根据内部状态获取 SOC 的测量值。 
%该函数假设所有状态都是可测量的。



% 检查是否使用了菲克扩散定律。 这是定义评估 SOC 的正确方法所必需的。
if(param{i}.SolidPhaseDiffusion~=3)
    cs_average = cs_average_t{i}(end,param{i}.Np+1:end);
else
    start_index = param{i}.Nr_p*param{i}.Np+1;
    end_index   = start_index+param{i}.Nr_n-1;
    cs_average  = zeros(param{i}.Nn,1);
    for n=1:param{i}.Nn
        cs_average(n)   = 1/param{i}.Rp_n*(param{i}.Rp_n/param{i}.Nr_n)*sum(cs_average_t{i}(end,start_index:end_index));
        start_index     = end_index + 1;
        end_index       = end_index + param{i}.Nr_n;
    end
end
Csout  = sum(cs_average);
Sout   = 100*((1/param{i}.len_n*(param{i}.len_n/(param{i}.Nn))*Csout/param{i}.cs_maxn)-param{i}.theta_min_neg)/(param{i}.theta_max_neg-param{i}.theta_min_neg);
end