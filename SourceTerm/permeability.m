%Define permeability%
function k_ = permeability(node)
%give permeability value at each node
    if node < 100
        k_ = 1.28e-14; %m^2
    else
        k_ = 1.28e-14; %m^2
    end
end