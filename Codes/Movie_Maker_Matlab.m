clc
clear

disp('------------------------------------------')
disp('|Beware, this code is for Matlab ONLY !!! |')
disp('------------------------------------------')

vidfile = VideoWriter('Output.mp4','MPEG-4');
open(vidfile);
listing = dir('*.png');
for i=1:1:length(listing)
    name=listing(i).name
    frame=imread(name);
    [height, width, null]=size(frame);
    frame=imresize(frame,800/height,'nearest');
    writeVideo(vidfile, frame);
end
close(vidfile)
