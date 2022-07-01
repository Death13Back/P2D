function [Y0_existence,YP0_existence,Y0,YP0] = checkInitialStates(initialState)
% checkInitialStates 检查用户提供的初始状态结构是否合适。



if(~isempty(initialState))
    Y_field_existence        = isfield(initialState,'Y');
    YP_field_existence       = isfield(initialState,'YP');

    if(Y_field_existence==0 || YP_field_existence==0)
        clc
        error('When defining initial states structure, please remember to set the Y and YP field')
    end

    if(Y_field_existence==1 && YP_field_existence==1)
        if(~isempty(initialState.Y) && ~isempty(initialState.Y))
            Y0              = initialState.Y;
            YP0             = initialState.YP;
            Y0_existence    = 1;
            YP0_existence   = 1;
        else
            Y0              = [];
            YP0             = [];
            Y0_existence    = 0;
            YP0_existence   = 0;
        end
    end
else
    Y0_existence    = 0;
    YP0_existence   = 0;
    Y0              = [];
    YP0             = [];
end
end
