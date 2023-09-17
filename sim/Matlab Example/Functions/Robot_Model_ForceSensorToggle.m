function [] = Robot_Model_ForceSensorToggle(iR,state)
global DTL

DTL.Robot{iR}.ForceSensor.State = state;

if state==1

        DTL.Robot{iR}.ForceSensor.X.Visible = 'off';
        DTL.Robot{iR}.ForceSensor.Y.Visible = 'off';
        DTL.Robot{iR}.ForceSensor.Z.Visible = 'off';
        DTL.Robot{iR}.ForceSensor.Mag.Visible = 'on';
        
else
    
        DTL.Robot{iR}.ForceSensor.X.Visible = 'off';
        DTL.Robot{iR}.ForceSensor.Y.Visible = 'off';
        DTL.Robot{iR}.ForceSensor.Z.Visible = 'off';
        DTL.Robot{iR}.ForceSensor.Mag.Visible = 'off';    
end

end