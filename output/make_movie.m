v = VideoWriter('newfile');
%v.CompressionRatio = 3;
%v.Quality = 80;
open(v)
for i = 1 : nTstep
    A = imread(sprintf('snapshots/step%04i.png',i));
    writeVideo(v,A)
end
close(v)