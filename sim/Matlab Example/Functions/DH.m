function [T] = DH(alp,a,d,the)

T = [cos(the)          -sin(the)           0           a;
    sin(the)*cos(alp)   cos(the)*cos(alp)   -sin(alp)   -sin(alp)*d;
    sin(the)*sin(alp)   cos(the)*sin(alp)   cos(alp)    cos(alp)*d;
    0                   0                   0           1];

end