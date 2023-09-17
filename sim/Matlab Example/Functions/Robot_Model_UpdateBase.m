function [] = Robot_Model_UpdateBase(iR,T,ax1)
global DTL
Config = DTL.Robot{iR}.Config;

%% DH

L1 = DTL.Robot{iR}.FixedLengths(1);  
L2 = DTL.Robot{iR}.FixedLengths(2);
L3 = DTL.Robot{iR}.FixedLengths(3);
L4 = DTL.Robot{iR}.FixedLengths(4);

th1 = Config(1);
th2 = Config(2);
th3 = Config(3);
th4 = Config(4);
th5 = Config(5);
th6 = Config(6);
th7 = Config(7);
DTL.Robot{iR}.TW_0=T;
DTL.Robot{iR}.T0_{8} = T;
DTL.Robot{iR}.T0_{1} = T*DH(0, 0, L1, deg2rad(th1));
DTL.Robot{iR}.T0_{2} = DTL.Robot{iR}.T0_{1}*DH(-pi/2, 0, 0, deg2rad(th2));
DTL.Robot{iR}.T0_{3} = DTL.Robot{iR}.T0_{2}*DH(pi/2, 0, L2, deg2rad(th3));
DTL.Robot{iR}.T0_{4} = DTL.Robot{iR}.T0_{3}*DH(pi/2, 0, 0, deg2rad(th4));
DTL.Robot{iR}.T0_{5} = DTL.Robot{iR}.T0_{4}*DH(-pi/2, 0, L3, deg2rad(th5));
DTL.Robot{iR}.T0_{6} = DTL.Robot{iR}.T0_{5}*DH(-pi/2, 0, 0, deg2rad(th6));
DTL.Robot{iR}.T0_{7} = DTL.Robot{iR}.T0_{6}*DH(pi/2, 0, L4, deg2rad(th7));

for iJ = 1:7
    DTL.Robot{iR}.Joint{iJ}.Transform.Matrix = DTL.Robot{iR}.T0_{iJ} * DTL.Robot{iR}.Joint{iJ}.Offset.Matrix;
end
DTL.Robot{iR}.Joint{8}.Transform.Matrix = DTL.Robot{iR}.T0_{8};
Robot_Model_updateAxes(iR,ax1)
end