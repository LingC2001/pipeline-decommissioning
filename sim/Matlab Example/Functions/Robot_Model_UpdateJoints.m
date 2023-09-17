function [] = Robot_Model_UpdateJoints(iR, th1, th2, th3, th4, th5, th6, th7, ax1)
global DTL

%% DH

L1 = DTL.Robot{iR}.FixedLengths(1);  
L2 = DTL.Robot{iR}.FixedLengths(2);
L3 = DTL.Robot{iR}.FixedLengths(3);
L4 = DTL.Robot{iR}.FixedLengths(4);


%Joint limit 
th1 = min(max((th1),(-170)),(170));
th2 = min(max((th2),(-120)),(120));
th3 = min(max((th3),(-170)),(170));
th4 = min(max((th4),(-120)),(120));
th5 = min(max((th5),(-170)),(170));
th6 = min(max((th6),(-120)),(120));
th7 = min(max((th7),(-175)),175);

DTL.Robot{iR}.Config = [th1, th2, th3, th4, th5, th6, th7];

DTL.Robot{iR}.T0_{8} = DTL.Robot{iR}.TW_0;
DTL.Robot{iR}.T0_{1} = DTL.Robot{iR}.TW_0*DH(0, 0, L1, deg2rad(th1));
DTL.Robot{iR}.T0_{2} = DTL.Robot{iR}.T0_{1}*DH(-pi/2, 0, 0, deg2rad(th2));
DTL.Robot{iR}.T0_{3} = DTL.Robot{iR}.T0_{2}*DH(pi/2, 0, L2, deg2rad(th3));
DTL.Robot{iR}.T0_{4} = DTL.Robot{iR}.T0_{3}*DH(pi/2, 0, 0, deg2rad(th4));
DTL.Robot{iR}.T0_{5} = DTL.Robot{iR}.T0_{4}*DH(-pi/2, 0, L3, deg2rad(th5));
DTL.Robot{iR}.T0_{6} = DTL.Robot{iR}.T0_{5}*DH(-pi/2, 0, 0, deg2rad(th6));
DTL.Robot{iR}.T0_{7} = DTL.Robot{iR}.T0_{6}*DH(pi/2, 0, L4, deg2rad(th7));

for iJ = 1:7
    set(DTL.Robot{iR}.Joint{iJ}.Transform,'Matrix',DTL.Robot{iR}.T0_{iJ} * DTL.Robot{iR}.Joint{iJ}.Offset.Matrix)
end
set(DTL.Robot{iR}.Joint{8}.Transform,'Matrix',DTL.Robot{iR}.T0_{8});

if DTL.Robot{iR}.AxesProp.state==1||DTL.Robot{iR}.AxesProp.state==-1
Robot_Model_updateAxes(iR,ax1)
end

if DTL.Robot{iR}.ForceSensor.State==1
Robot_Model_updateForceSensor(iR,ax1)
end

end