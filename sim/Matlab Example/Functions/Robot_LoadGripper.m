function [] = Robot_LoadGripper(iR,STL,ColourF,ColourE,AlphaF,AlphaE,WidthE,compressed,ax1)
global DTL
% MAKE SURE THE ORIGNAL STL HAS THE CORRECT ORIGIN 
% Attached to Gripper Frame

if compressed==1
A = stlread(append("STL/gripper/compressed/simplify_",STL,".stl"));
else
A = stlread(append("STL/gripper/",STL,".stl")); 
end
R = [0 -1 0; 1 0 0; 0 0 1];
B = R*A.Points';

DTL.Robot{iR}.Gripper.stl = triangulation(A.ConnectivityList, B');

DTL.Robot{iR}.Gripper.Transform = hgtransform(ax1); %gg
DTL.Robot{iR}.Gripper.Offset = hgtransform(ax1);
DTL.Robot{iR}.Gripper.fig=trimesh(DTL.Robot{iR}.Gripper.stl,'Parent',ax1,'Parent',DTL.Robot{iR}.Joint{7}.Transform); %%hh
DTL.Robot{iR}.Gripper.fig.FaceColor = ColourF;
DTL.Robot{iR}.Gripper.fig.EdgeColor = ColourE;
DTL.Robot{iR}.Gripper.fig.FaceAlpha = AlphaF;
DTL.Robot{iR}.Gripper.fig.EdgeAlpha = AlphaE;
DTL.Robot{iR}.Gripper.fig.LineWidth = WidthE;
DTL.Robot{iR}.Gripper.fig.FaceLighting = "flat";

end