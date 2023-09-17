function [] = Robot_Model_updateForceSensor(iR,ax1)
global DTL

%DTL.Robot{iR}.ForceSensor.F = DTL.Robot{iR}.ForceSensor.F -0.5+[rand,rand,rand]*1;
if DTL.Robot{iR}.ForceSensor.State==1

    Length = DTL.Robot{iR}.ForceSensor.Length;
    F = DTL.Robot{iR}.ForceSensor.F;
    
    R_base = DTL.Robot{iR}.TW_0(1:3,1:3);

    F = R_base*F';

    T = DTL.Robot{iR}.T0_{7};
    p = T(1:3,4);
    Fs = F*Length;

    set(DTL.Robot{iR}.ForceSensor.X,'XData',[p(1) p(1)+Fs(1)],'YData',[p(2) p(2)],'ZData',[p(3) p(3)]);
    set(DTL.Robot{iR}.ForceSensor.Y,'XData',[p(1) p(1)],'YData',[p(2) p(2)+Fs(2)],'ZData',[p(3) p(3)]);
    set(DTL.Robot{iR}.ForceSensor.Z,'XData',[p(1) p(1)],'YData',[p(2) p(2)],'ZData',[p(3),p(3)+Fs(3)]);
    set(DTL.Robot{iR}.ForceSensor.Mag,'XData',[p(1) p(1)+Fs(1)],'YData',[p(2) p(2)+Fs(2)],'ZData',[p(3) p(3)+Fs(3)]);

end

end