function [] = Lab_LoadRobot(iR,ColourF,ColourE,AlphaF,AlphaE,WidthE,model,ax1)
global DTL

if (model == 0) || (model == 14)
    L1 = 0.360;
    L2 = 0.420;
    L3 = 0.400;
    L4 = 0.152-0.0006;
    DTL.Robot{iR}.FixedLengths = [L1 L2 L3 L4];

    DTL.Robot{iR}.TW_0 = double(eye(4));
    DTL.Robot{iR}.T0_{1} = DTL.Robot{iR}.TW_0*DH(0, 0, L1, 0);
    DTL.Robot{iR}.T0_{2} = DTL.Robot{iR}.T0_{1}*DH(-pi/2, 0, 0, 0);
    DTL.Robot{iR}.T0_{3} = DTL.Robot{iR}.T0_{2}*DH(pi/2, 0, L2, 0);
    DTL.Robot{iR}.T0_{4} = DTL.Robot{iR}.T0_{3}*DH(pi/2, 0, 0, 0);
    DTL.Robot{iR}.T0_{5} = DTL.Robot{iR}.T0_{4}*DH(-pi/2, 0, L3, 0);
    DTL.Robot{iR}.T0_{6} = DTL.Robot{iR}.T0_{5}*DH(-pi/2, 0, 0, 0);
    DTL.Robot{iR}.T0_{7} = DTL.Robot{iR}.T0_{6}*DH(pi/2, 0, L4, 0);
    DTL.Robot{iR}.T0_{8} = DTL.Robot{iR}.TW_0;

    DTL.Robot{iR}.Config = zeros(7,1);

    offset1 = [0 0 0.1575-L1 0 0 0];
    offset2 = [0 0 0 0 0 pi];
    offset3 = [0 0 0.2045-L2 0 0 0];
    offset4 = [0 0 0 0 0 0];
    offset5 = [0 0 0.1845-L3 0 0 pi];
    offset6 = [0 0 0 0 0 pi];
    offset7 = [0 0 0 0 0 0];

elseif (model == 7)
    L1 = 0.340;
    L2 = 0.400;
    L3 = 0.400;
    L4 = 0.152-0.0006;
    DTL.Robot{iR}.FixedLengths = [L1 L2 L3 L4];


    DTL.Robot{iR}.TW_0 = double(eye(4));
    DTL.Robot{iR}.T0_{1} = DTL.Robot{iR}.TW_0*DH(0, 0, L1, 0);
    DTL.Robot{iR}.T0_{2} = DTL.Robot{iR}.T0_{1}*DH(-pi/2, 0, 0, 0);
    DTL.Robot{iR}.T0_{3} = DTL.Robot{iR}.T0_{2}*DH(pi/2, 0, L2, 0);
    DTL.Robot{iR}.T0_{4} = DTL.Robot{iR}.T0_{3}*DH(pi/2, 0, 0, 0);
    DTL.Robot{iR}.T0_{5} = DTL.Robot{iR}.T0_{4}*DH(-pi/2, 0, L3, 0);
    DTL.Robot{iR}.T0_{6} = DTL.Robot{iR}.T0_{5}*DH(-pi/2, 0, 0, 0);
    DTL.Robot{iR}.T0_{7} = DTL.Robot{iR}.T0_{6}*DH(pi/2, 0, L4, 0);
    DTL.Robot{iR}.T0_{8} = DTL.Robot{iR}.TW_0;

    DTL.Robot{iR}.Config = zeros(7,1);

    offset1 = [0 0 0.1575-L1 0 0 0];
    offset2 = [0 0 0 0 0 pi];
    offset3 = [0 0 0.1845-L2 0 0 0];
    offset4 = [0 0 0 0 0 0];
    offset5 = [0 0 0.1845-L3 0 0 pi];
    offset6 = [0 0 0 0 0 pi];
    offset7 = [0 0 0 0 0 0];
end

Robot_Model_LoadJoint(iR,8,ColourF,ColourE,AlphaF,AlphaE,WidthE,ax1,[0 0 0 0 0 0],model)
Robot_Model_LoadJoint(iR,1,ColourF,ColourE,AlphaF,AlphaE,WidthE,ax1,offset1,model)
Robot_Model_LoadJoint(iR,2,ColourF,ColourE,AlphaF,AlphaE,WidthE,ax1,offset2,model)
Robot_Model_LoadJoint(iR,3,ColourF,ColourE,AlphaF,AlphaE,WidthE,ax1,offset3,model)
Robot_Model_LoadJoint(iR,4,ColourF,ColourE,AlphaF,AlphaE,WidthE,ax1,offset4,model)
Robot_Model_LoadJoint(iR,5,ColourF,ColourE,AlphaF,AlphaE,WidthE,ax1,offset5,model)
Robot_Model_LoadJoint(iR,6,ColourF,ColourE,AlphaF,AlphaE,WidthE,ax1,offset6,model)
Robot_Model_LoadJoint(iR,7,ColourF,ColourE,AlphaF,AlphaE,WidthE,ax1,offset7,model)

Robot_Model_LoadAxes(iR,[1:8],0.15,2,ax1)
Robot_Model_AxesToggle(iR,[1:8],0,-1)

end