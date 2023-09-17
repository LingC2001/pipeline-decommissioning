function [] = Robot_Model_LoadJoint(iR,iJ,ColourF,ColourE,AlphaF,AlphaE,WidthE,ax1,offset,model)
global DTL

if iJ==8
    iJtext = '0';
else
    iJtext = string(iJ);
end

%Import STL
if model==7
DTL.Robot{iR}.Joint{iJ}.stl = stlread(append("STL/iiwa7/visual/link_",iJtext,".stl")); %8=B for Base
else
DTL.Robot{iR}.Joint{iJ}.stl = stlread(append("STL/iiwa14/visual/link_",iJtext,".stl")); %8=B for Base
end
%Plot STL
DTL.Robot{iR}.Joint{iJ}.Transform = hgtransform(ax1);
DTL.Robot{iR}.Joint{iJ}.Offset = hgtransform(ax1);
DTL.Robot{iR}.Joint{iJ}.fig=trimesh(DTL.Robot{iR}.Joint{iJ}.stl,'Parent',ax1,'Parent',DTL.Robot{iR}.Joint{iJ}.Transform);
DTL.Robot{iR}.Joint{iJ}.fig.FaceColor = ColourF;
DTL.Robot{iR}.Joint{iJ}.fig.EdgeColor = ColourE;
DTL.Robot{iR}.Joint{iJ}.fig.FaceAlpha = AlphaF;
DTL.Robot{iR}.Joint{iJ}.fig.EdgeAlpha = AlphaE;
DTL.Robot{iR}.Joint{iJ}.fig.LineWidth = WidthE;
DTL.Robot{iR}.Joint{iJ}.fig.FaceLighting = "flat";





DTL.Robot{iR}.Joint{iJ}.Transform.Matrix = DTL.Robot{iR}.Joint{iJ}.Transform.Matrix * makehgtform('translate',offset(1:3));
DTL.Robot{iR}.Joint{iJ}.Transform.Matrix = DTL.Robot{iR}.Joint{iJ}.Transform.Matrix * makehgtform('xrotate',offset(4));
DTL.Robot{iR}.Joint{iJ}.Transform.Matrix = DTL.Robot{iR}.Joint{iJ}.Transform.Matrix * makehgtform('yrotate',offset(5));
DTL.Robot{iR}.Joint{iJ}.Transform.Matrix = DTL.Robot{iR}.Joint{iJ}.Transform.Matrix * makehgtform('zrotate',offset(6));
DTL.Robot{iR}.Joint{iJ}.Offset.Matrix = DTL.Robot{iR}.Joint{iJ}.Transform.Matrix;
DTL.Robot{iR}.Joint{iJ}.Transform.Matrix = DTL.Robot{iR}.T0_{iJ} * DTL.Robot{iR}.Joint{iJ}.Offset.Matrix;
DTL.Robot{iR}.Home = DTL.Robot{iR}.Joint{iJ}.Transform.Matrix;



end