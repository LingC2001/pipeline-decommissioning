function [] = Robot_Model_LoadAxes(iR,selector,Length,WidthL,ax1)
global DTL

DTL.Robot{iR}.selector.loaded = selector;
DTL.Robot{iR}.selector.visible = [];
DTL.Robot{iR}.AxesProp.Length = Length;
DTL.Robot{iR}.AxesProp.Width = WidthL;

for i = selector

    T = DTL.Robot{iR}.T0_{i};
    p = T(1:3,4);
    X = T(1:3,1)*Length;
    Y = T(1:3,2)*Length;
    Z = T(1:3,3)*Length;

    DTL.Robot{iR}.Joint{i}.axes.X = plot3([p(1) p(1)+X(1)],[p(2) p(2)+X(2)],[p(3) p(3)+X(3)],'-r','Parent',ax1);
    DTL.Robot{iR}.Joint{i}.axes.Y = plot3([p(1) p(1)+Y(1)],[p(2) p(2)+Y(2)],[p(3) p(3)+Y(3)],'-g','Parent',ax1);
    DTL.Robot{iR}.Joint{i}.axes.Z = plot3([p(1) p(1)+Z(1)],[p(2) p(2)+Z(2)],[p(3) p(3)+Z(3)],'-b','Parent',ax1);
    DTL.Robot{iR}.Joint{i}.axes.X.LineWidth = WidthL;
    DTL.Robot{iR}.Joint{i}.axes.Y.LineWidth = WidthL;
    DTL.Robot{iR}.Joint{i}.axes.Z.LineWidth = WidthL;

end


end