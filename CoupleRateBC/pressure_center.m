function [p_input] = pressure_center(t,t_start,t_end)
%Define fluid injection (pressure change)
%Instead of using heaviside function, define tanh function
max_p = 2e6;
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
        p_input = 1e6 * abs(tanh((t-410))+1);
    else
        p_input = 1e6 * abs(tanh(-(t-710))+1);
    end

end
