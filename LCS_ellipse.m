function []  = LCS_ellipse()
%% detect ellipses in stereo images
%%
% read images
load('input.mat',...
    'sphImagePath1', 'sphImagePath2', 'sphImageNames');
[imgs1, imgs2] = readStereoImgs(sphImagePath1, sphImageNames, sphImagePath2, sphImageNames);
numImgs = length(imgs1);
ellipses1 = zeros(numImgs, 5);
ellipses2 = zeros(numImgs, 5);

%%
%detect ellipses
for i = 1 : numImgs
    disp('--------------------------------------------------------------');
    disp(i);
    disp('--------------------------------------------------------------');
    ellipseTmp = LCS_ellipse_from_img(imgs1{i});
    ellipses1(i, :) = ellipseTmp(1, :);
    ellipseTmp = LCS_ellipse_from_img(imgs2{i});
    ellipses2(i, :) = ellipseTmp(1, :);
end

%%
% save result
save('output.mat', 'ellipses1', 'ellipses2');
end

function [ellipses]  = LCS_ellipse_from_img(I)
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

% parameters
Tac = 165;
Tr = 0.6;
specified_polarity = 0;

%%
% read image 
%image path
%filename = 'C:\Users\totht\ELTE\SphereCalib\examples\example4\Scene2\Sphere\2D1\(1).jpg'
%'C:\Users\totht\ELTE\High-quality-ellipse-detection\pics\666.jpg';
%disp('------read image------');
%I = imread(filename);


%% detecting ellipses from real-world images
[ellipses, ~, posi] = ellipseDetectionByArcSupportLSs(I, Tac, Tr, specified_polarity);

disp('draw detected ellipses');
drawEllipses(ellipses',I);
% display
ellipses(:,5) = ellipses(:,5)./pi*180;
disp(ellipses);
disp(['The total number of detected ellipses£º',num2str(size(ellipses,1))]);

%% draw ellipse centers
%hold on;
%candidates_xy = round(posi+0.5);%candidates' centers (col_i, row_i)
%plot(candidates_xy(:,1),candidates_xy(:,2),'.');%draw candidates' centers.

%% write the result image
%set(gcf,'position',[0 0 size(I,2) size(I,1)]);
%saveas(gcf, 'D:\Graduate Design\Ellipse Detection\MyEllipse - github\pics\666_all.jpg', 'jpg');
end

function [imgs1, imgs2] = readStereoImgs(imgDir1, imgNames1, imgDir2, imgNames2)
%
% Read stereo image
%
% 
% INPUT
% Source file infos
% - imgDir1, imgNames1 -- first image path and file names
% - imgDir2, imgNames2 -- second image path and file names
%
%
numImgs1 = length(imgNames1);
if numImgs1 < 1
    error('Not enough input image file names (second argument).');
end

numImgs2 = length(imgNames2);
if numImgs2 < 1
    error('Not enough input image file names (fourth argument).');
end

if numImgs1 ~= numImgs2
    error('Number of input stereo images are different (second and fourth argument).');
end

imgs1 = readImgs(imgDir1, imgNames1);
imgs2 = readImgs(imgDir2, imgNames2);

end

function [imgs] = readImgs(imgDir, imgNames)
%
% Read image
%
% 
% INPUT
% Source file infos
% - imgDir, imgNames -- image path and file names
%
%
imgs = cell(length(imgNames), 1);

for i= 1:length(imgNames)
    imgFile = strcat(imgDir,imgNames(i));
    imgs{i} = imread(char(imgFile));
end

end



