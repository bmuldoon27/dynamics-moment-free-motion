currentDir = cd;
picDir = [currentDir '/gif_figures_asym'];
cd(picDir);
gifname = 'mfmrb_asym.gif';
imind = zeros(1006, 1382, 1, 9);

f1 = figure(1);
for iii = 1:400
    filename = ['iter=' num2str(iii) '.000000.png'];
    imshow(filename);
    drawnow;
    
    frame = getframe(f1);
    im = frame2im(frame);
    [imind,cm]=rgb2ind(im,256);
    if iii == 1
        imwrite(imind,cm,gifname,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,gifname,'gif','WriteMode','append');
    end
end

cd(currentDir);