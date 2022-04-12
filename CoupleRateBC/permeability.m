%Define permeability%
function k_ = permeability(node)
%give permeability value at each node
    if abs(node) < 100
        k_ = 1.28e-14; %m^2
    else
        k_ = 1.28e-14; %m^2
%         k_ = 1.28e-14;
    end
end