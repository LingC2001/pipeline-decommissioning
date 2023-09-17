clear all
close all
addpath('Functions/')
global DTL


%% Variables
ColourF=[0.4660 1 0.1880];
%ColourE=[0.4660 0.5 0.1880];
ColourE='none';
WidthE=0.4;
AlphaF=0.25;
AlphaE=0.8;
axRGBl = 0.2;

f1 = figure(1); 
ax1 = axes(f1); 
f1.WindowState = 'maximized';

hold(ax1,'on') 
axis(ax1,'equal') 
light(ax1,'Position',[100 100 100],'Style','local') 
grid(ax1,'on')
view(ax1,30,15)
% xlim([-1.8 1.8])
% ylim([-0.8 0.8])
% zlim([-0.2 1])

%% Load Robots
Lab_LoadRobot(1,'red','red',0.2,AlphaE,WidthE,7,ax1)
%Lab_LoadRobot(3,'green','green',0.1,0.1,WidthE,0,ax1)
%Lab_LoadRobot(4,'blue',ColourE,AlphaF,AlphaE,WidthE,0,ax1)

%% Load Grippers
Robot_LoadGripper(1,'2F85','red',ColourE,1,AlphaE,WidthE,0,ax1)
%Robot_LoadGripper(3,'2F140','black',ColourE,AlphaF,AlphaE,WidthE,0,ax1)
%Robot_LoadGripper(4,'2F85','black',ColourE,AlphaF,AlphaE,WidthE,0,ax1)

%% Load Force Sensor 
Robot_Model_LoadForceSensor(1,0.1,5,[0,0,0],ax1)
%Robot_Model_LoadForceSensor(3,0.4,3,[0,0,0],ax1)
%Robot_Model_LoadForceSensor(4,0.1,3,[0,0,0],ax1)
 
Robot_Model_ForceSensorToggle(1,1)
%Robot_Model_ForceSensorToggle(3,0)
%Robot_Model_ForceSensorToggle(4,0)

%% Update Base position
%T3 = [-1 0 0 -1.35; 0 -1 0 -0.015; 0 0 1 0.1; 0 0 0 1];
%T4 = [-1 0 0 1.7; 0 -1 0 -0.3; 0 0 1 0.1; 0 0 0 1];
Robot_Model_UpdateBase(1,[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1], ax1)
%Robot_Model_UpdateBase(3,T3, ax1)
%Robot_Model_UpdateBase(4,T4, ax1)

%TO DO: Add in base anaimations

%% Toggle Bases and EE Axes
Robot_Model_AxesToggle(1,[7 8],1,-1,ax1)
Robot_Model_AxesToggle(1,[7 8],1,-1,ax1)
%Robot_Model_AxesToggle(3,[7 8],1,1,ax1)
%Robot_Model_AxesToggle(4,[7 8],0,0,ax1)
 
%% Build Lab
% TO DO: turn shelving object to smart storage function 
% Object 1 Main Benchtop
T_o1 = [1 0 0 0;
    0 1 0 0.180;
    0 0 1 0;
    0 0 0 1];
Lab_LoadObject(1,0,T_o1,'Benchtop_Square','#858585','none',0.8,1,1,ax1)
%{
% Object 3 FlexFellow
T_o3 = DTL.Robot{3}.TW_0;
T_o3(1:3,1:3) = [-1 0 0; 0 -1 0; 0 0 1]; 
T_o3(3,4) = -0.85;
Lab_LoadObject(3,0,T_o3,'FlexFellow','#858585','none',0.8,1,1,ax1)

% Object 4 FlexFellow
T_o4 = DTL.Robot{4}.TW_0;
T_o4(1:3,1:3) = [-1 0 0; 0 -1 0; 0 0 1]; 
T_o4(3,4) = -0.85;
Lab_LoadObject(4,0,T_o4,'FlexFellow','#858585','none',0.8,1,1,ax1)

% Object 2 TV 
Lab_LoadObject(4,1,[-1 0.9 0.05],[2 0.01 1.1],'black','black',0.2,1,WidthE,ax1) %FF2

% Object 5 Shelving 
T_o5 = [0 0 -1 -2.8; -1 0 0 1.3; 0 1 0 -0.8; 0 0 0 1];
Lab_LoadObject(5,0,T_o5,'Shelving1','#858585','none',0.8,1,1,ax1)

% Object 6 Rectangle Benchtop
T_o6 = [0 1 0 2;
    -1 0 0 0.7;
    0 0 1 0;
    0 0 0 1];
Lab_LoadObject(6,0,T_o6,'Benchtop_Rectangle','#858585','none',0.8,1,1,ax1)

% Object 7 Printer 3
T_o7a = [-1 0 0 0;
    0 -1 0 0;
    0 0 1 0;
    0 0 0 1];
T_o7b = [1 0 0 -0.75;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];
Lab_LoadObject(7,0,T_o7b*T_o6*T_o7a,'Ender3-V2','yellow','none',0.8,1,1,ax1)

% Object 8 Printer 4
T_o8a = [-1 0 0 0;
    0 -1 0 0;
    0 0 1 0;
    0 0 0 1];
T_o8b = [1 0 0 -0.2;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];
Lab_LoadObject(8,0,T_o8b*T_o6*T_o8a,'Ender3-V2','yellow','none',0.8,1,1,ax1)

% Walls
Lab_LoadObject(9,1,[-2.8 -2 -0.801],[6.7 3.2 0.01],[0.4660 0.6740 0.1880],'black',0.1,0.4,0.2,ax1) %Floor
Lab_LoadObject(10,1,[-2.8 1.3 -0.80],[6.7 0.01 3],[0.4660 0.6740 0.1880],'black',0.1,0.4,0.2,ax1) %Wall Back
Lab_LoadObject(11,1,[3.9 -2 -0.80],[0.01 3.2 3],[0.4660 0.6740 0.1880],'black',0.1,0.4,0.2,ax1) %Wall Side
% 

% Object 12 Rectangle Benchtop
T_o12 = [1 0 0 3.5;
    0 1 0 -0.6;
    0 0 1 0;
    0 0 0 1];
Lab_LoadObject(12,0,T_o12,'Benchtop_Rectangle','#858585','none',0.8,1,1,ax1)

% Objects 13+
AA = load('FilamentFake.mat');
c = length(AA.FilamentFake.Color);

for i =(1:c)
    ii = i+12;
    T = [1 0 0 0.6; 0 1 0 -2.2+((i-1)*0.15); 0 0 1 2.05; 0 0 0 1]*T_o5;
    Lab_LoadObject(ii,0,T,'FilamentSpool',AA.FilamentFake.Color{i},'none',0.8,1,0.4,ax1);
end


% Objects 27+
BB = load('PrintedPartsFake.mat');
c = length(BB.PrintedPartsFake.Color);
for i = 1:c
    ii = i+26;
    x = BB.PrintedPartsFake.LocationX(i);
    y = BB.PrintedPartsFake.LocationY(i);
    level = BB.PrintedPartsFake.Level(i);
    T = T_shelf(level,x,y)*T_o5*[1 0 0 0; 0 0 1 0; 0 -1 0 0; 0 0 0 1]*[-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1];
    name = BB.PrintedPartsFake.PartName{i};
    colour = BB.PrintedPartsFake.Color{i};
    if strcmp(colour,'purple')
        colour = [0.5 0 0.5];
    elseif strcmp(colour,'brown')
        colour = [0.588 0.294 0];
    end
    Lab_LoadObject(ii,0,T,name,colour,'none',0.8,1,0.4,ax1);
end
%}


%% Animate

NewPose1 = [-82, 16, -9, -114, -5, -42, 3];
%NewPose2 = [9, -41, -31, 56, 22, -89, 82];
%NewPose3 = [-62, 25, 63, 117, 20, 97, -83];

%NewPose = [NewPose1;NewPose2;NewPose3];
res = 100;
%Check and update function
 %Robot_Model_JSPanimate([1 3 4],NewPose,res,ax1)
Robot_Model_JSPanimate([1],NewPose1,res,ax1)