function [] = Robot_Model_updateAxes(iR,ax1)
global DTL

selectorV = DTL.Robot{iR}.selector.visible;
Length = DTL.Robot{iR}.AxesProp.Length;

for i = selectorV

    T = DTL.Robot{iR}.T0_{i};
    p = T(1:3,4);
    X = T(1:3,1)*Length;
    Y = T(1:3,2)*Length;
    Z = T(1:3,3)*Length;

    set(DTL.Robot{iR}.Joint{i}.axes.X,'XData',[p(1) p(1)+X(1)],'YData',[p(2) p(2)+X(2)],'ZData',[p(3) p(3)+X(3)]);
    set(DTL.Robot{iR}.Joint{i}.axes.Y,'XData',[p(1) p(1)+Y(1)],'YData',[p(2) p(2)+Y(2)],'ZData',[p(3) p(3)+Y(3)]);
    set(DTL.Robot{iR}.Joint{i}.axes.Z,'XData',[p(1) p(1)+Z(1)],'YData',[p(2) p(2)+Z(2)],'ZData',[p(3) p(3)+Z(3)]);
    
end

end