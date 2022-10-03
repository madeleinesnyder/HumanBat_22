cap = VideoReader('head_dir_tracking/head_dir_tracking.mp4');
frame = read(cap, 500);
tiledlayout(2,4);
nexttile;
imshow(frame)
frame0 = read(cap,1);
nexttile;
imshow(frame0)
nexttile;
diff = frame-frame0;
imshow(diff)
nexttile;
mask = abs(diff) > 120;
masked = immultiply(diff,mask);
imshow(masked)
nexttile;
bin = imbinarize(rgb2gray(masked), 0.4);
imshow(bin)
BW = imclearborder(bin);
BW = bwareafilt(BW,1);
BW = bwconvhull(BW);
nexttile;
imshow(BW)
stats = regionprops(BW,'Centroid','MajorAxisLength','MinorAxisLength', 'Orientation')
drawellipse('Center', stats.Centroid, 'RotationAngle', stats.Orientation, 'SemiAxes', [stats.MajorAxisLength/2 stats.MinorAxisLength/2], 'MarkerSize',0.1)
nexttile;
imshow(frame)
drawellipse('Center', stats.Centroid, 'RotationAngle', stats.Orientation, 'SemiAxes', [stats.MajorAxisLength/2 stats.MinorAxisLength/2], 'MarkerSize',0.1)
