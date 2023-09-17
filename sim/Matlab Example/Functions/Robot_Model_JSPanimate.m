function [] = Robot_Model_JSPanimate(iR,NewPose,res,ax1)

global DTL
% iR should be a 1xn vector indicating target robots
% New Pose should be a nx7 Matrix

L = length(iR);

CurPose = zeros(L,7);
JSP = zeros(res,7,L);

%% Load Current Pose

for i=1:L
    CurPose(i,:)=DTL.Robot{iR(i)}.Config;
    JSP(:,:,i) = [linspace(CurPose(i,1),NewPose(i,1),res);
        linspace(CurPose(i,2),NewPose(i,2),res);
        linspace(CurPose(i,3),NewPose(i,3),res);
        linspace(CurPose(i,4),NewPose(i,4),res);
        linspace(CurPose(i,5),NewPose(i,5),res);
        linspace(CurPose(i,6),NewPose(i,6),res);
        linspace(CurPose(i,7),NewPose(i,7),res);]';
end

for j=1:res
    for i=1:L
        Robot_Model_UpdateJoints(iR(i), JSP(j,1,i), JSP(j,2,i), JSP(j,3,i), JSP(j,4,i), JSP(j,5,i), JSP(j,6,i), JSP(j,7,i), ax1)
    end
    drawnow;
end

end