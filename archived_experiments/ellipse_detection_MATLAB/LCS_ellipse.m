function []  = LCS_ellipse()
%% parameters illustration
%1) Tac: 
%The threshold of elliptic angular coverage which ranges from 0~360. 
%The higher Tac, the more complete the detected ellipse should be.
%2) Tr:
%The ratio of support inliers to ellipse which ranges from 0~1.
%The higher Tr, the more sufficient the support inliers are.
%3) specified_polarity: 
%1 means detecting the ellipses with positive polarity;
%-1 means detecting the ellipses with negative polarity; 
%0 means detecting all ellipses from image


close all;

%image path
filename = 'C:\Users\asus\pipe_repo\pipe-recognisation\High-quality-ellipse-detection-master\pics\p23.jpg';
%filename = 'C:\Users\asus\pipe_repo\pipe-recognisation\High-quality-ellipse-detection-master\canny_images\c3.jpg';

% parameters
Tac = 165; %default 165
Tr = 0.6;  %default 0.6
specified_polarity = 0;
num_image_orientations = 1;

%%
% read image 
disp('------read image------');
img = imread(filename);
img = imresize(img, [1024,1024]);

%% detecting ellipses from real-world images
[ellipses, ~, posi] = ellipseDetectionByArcSupportLSs(img, Tac, Tr, specified_polarity);

disp('draw detected ellipses');
drawEllipses(ellipses',img);
% display
ellipses(:,5) = ellipses(:,5)./pi*180;
ellipses;
disp(['The total number of detected ellipses: ',num2str(size(ellipses,1))]);


% %% detecting ellipses from real-world images
% max_num_ellipses = -1;
% for i = 1:num_image_orientations
%     img_rot = imrotate(img, round((i-1)*360/num_image_orientations));
%     [ellipses, ~, posi] = ellipseDetectionByArcSupportLSs(img, Tac, Tr, specified_polarity);
% 
%     % display
%     ellipses(:,5) = ellipses(:,5)./pi*180;
%     ellipses
%     num_ellipses = size(ellipses,1);
%     if num_ellipses > max_num_ellipses
%         max_num_ellipses = num_ellipses;
%         max_ellipses = ellipses
%         max_posi = posi;
%         max_rot = img_rot;
%     end
%     disp(['The total number of detected ellipses: ',num2str(num_ellipses)]);
%     disp('-----------------------------------------------------------');
%     disp('-----------------------------------------------------------');
% end
% disp('draw detected ellipses');
% disp(['The final number of detected ellipses: ',num2str(max_num_ellipses)]);
% 
% drawEllipses(max_ellipses',img);

%% draw ellipse centers
%hold on;
%candidates_xy = round(max_posi+0.5);%candidates' centers (col_i, row_i)
%plot(candidates_xy(:,1),candidates_xy(:,2),'.');%draw candidates' centers.

%% write the result image
%set(gcf,'position',[0 0 size(I,2) size(max_rot,1)]);
%saveas(gcf, 'D:\Graduate Design\Ellipse Detection\MyEllipse - github\pics\666_all.jpg', 'jpg');
end



