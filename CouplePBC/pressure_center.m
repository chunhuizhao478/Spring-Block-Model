function [p_input,Dp_input] = pressure_center(t,t_start,t_end)
%Define fluid injection (pressure change)
%Instead of using heaviside function, define tanh function
%max_p duration: 500 - 1000
global max_p %fraction of normal stress 
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

        Dp_input = 0;
    elseif t_inject_start < t && t <= t_inject_end 
        p_input = max_p / 2 * abs(tanh(0.05*(t-450))+1);

        Dp_input = 0.025 * max_p * sech(22.5-0.05*t)^2;
    else
        p_input = max_p / 2 * abs(tanh(-0.05*(t-1050))+1);
        
        Dp_input = -0.025 * max_p * sech(0.05*(-1050 + t))^2;
    end

end
