%Define permeability%
function k_ = k(node)
%give permeability value at each node
    if abs(node) < 100
        k_ = 1.28e-14; %m^2
    else
        k_ = 1.28e-18; %m^2
    end
end