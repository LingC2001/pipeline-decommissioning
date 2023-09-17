function [] = Robot_Model_LoadForceSensor(iR,Length,WidthL,F,ax1)
global DTL

DTL.Robot{iR}.ForceSensor.Length = Length;
DTL.Robot{iR}.ForceSensor.Width = WidthL;
DTL.Robot{iR}.ForceSensor.State = 1;
DTL.Robot{iR}.ForceSensor.F = F;

T = DTL.Robot{iR}.T0_{7};
p = T(1:3,4);
Fs = F*Length;

DTL.Robot{iR}.ForceSensor.X = plot3([p(1) p(1)+Fs(1)],[p(2) p(2)],[p(3) p(3)],'-r','Parent',ax1);
DTL.Robot{iR}.ForceSensor.Y = plot3([p(1) p(1)],[p(2) p(2)+Fs(2)],[p(3) p(3)],'-g','Parent',ax1);
DTL.Robot{iR}.ForceSensor.Z = plot3([p(1) p(1)],[p(2) p(2)],[p(3) p(3)+Fs(3)],'-b','Parent',ax1);
DTL.Robot{iR}.ForceSensor.Mag = plot3([p(1) p(1)+Fs(1)],[p(2) p(2)+Fs(2)],[p(3) p(3)+Fs(3)],'-k','Parent',ax1);
DTL.Robot{iR}.ForceSensor.X.LineWidth = WidthL;
DTL.Robot{iR}.ForceSensor.Y.LineWidth = WidthL;
DTL.Robot{iR}.ForceSensor.Z.LineWidth = WidthL;
DTL.Robot{iR}.ForceSensor.Mag.LineWidth = WidthL;

end