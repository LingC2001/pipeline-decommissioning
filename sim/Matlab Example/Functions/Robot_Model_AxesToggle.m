function [] = Robot_Model_AxesToggle(iR,selector,state,stateall,ax1)
global DTL

if stateall == 1
    selector = 1:8;
    state = 1;
elseif stateall == 0
    selector = 1:8;
    state = 0;
else
    stateall = -1;
end

DTL.Robot{iR}.AxesProp.state = stateall;  

if state==1
    DTL.Robot{iR}.selector.visible = selector;
    Robot_Model_updateAxes(iR,ax1)
    for i = selector
        DTL.Robot{iR}.Joint{i}.axes.X.Visible = 'on';
        DTL.Robot{iR}.Joint{i}.axes.Y.Visible = 'on';
        DTL.Robot{iR}.Joint{i}.axes.Z.Visible = 'on';
        
    end
else
    for i = selector
        DTL.Robot{iR}.Joint{i}.axes.X.Visible = 'off';
        DTL.Robot{iR}.Joint{i}.axes.Y.Visible = 'off';
        DTL.Robot{iR}.Joint{i}.axes.Z.Visible = 'off';
    end
end

end