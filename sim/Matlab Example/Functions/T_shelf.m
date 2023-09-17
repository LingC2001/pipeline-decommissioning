function [T] = T_shelf(level,x,y)

T = eye(4);

x=x-1;
y=y-1;
level=level-1;

tray_w = 0.265;
tray_l = 0.538;

grid_x = 10;
grid_y = 20;

%% X Pos Calc
tray_x = tray_w/grid_x;
tray_gap = 0.028;
x_gap = floor(x/grid_x);
part_pos_x = x*tray_x + x_gap*tray_gap;

%% Y Pos Calc
tray_y = tray_l/grid_y;
part_pos_y = y*tray_y;

%% Z Pos Calc
part_pos_z = level*0.5;

%% 
T(3,4) = 0.66565 + 0.5 - part_pos_z;
T(2,4) = -0.091 - 2.318 + part_pos_x;
T(1,4) = 0.031 + 0.538 - part_pos_y;
end