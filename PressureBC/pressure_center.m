function [p_input] = pressure_center(t,t_start,t_end)
%Define fluid injection (pressure change)
%Instead of using heaviside function, define tanh function
global max_p
t_inject_start = t_start; 
t_inject_end = t_end;

% %option1
%     if t < t_inject_start
%         p_input = 0;
%     elseif t == t_inject_start
%         p_input = 1e6;
%     elseif t_inject_start < t && t < t_inject_end 
%         p_input = 2e6;
%     elseif t == t_inject_end
%         p_input = 1e6;
%     else
%         p_input = 0;
%     end

%option2
    if t <= t_inject_start
        p_input = 0;
    elseif t_inject_start < t && t <= t_inject_end 
%         p_input = max_p * abs(tanh(1e-2*(t-800))+1); 
        p_input = max_p / 2 * abs(tanh(0.05*(t-450))+1);
        %p_input = max_p; %%step func
    else
%         p_input = max_p * abs(tanh(-1e-2*(t-1900))+1);
        p_input = max_p / 2 * abs(tanh(-0.05*(t-1050))+1);
        %p_input = 0; %%step func
    end

end
