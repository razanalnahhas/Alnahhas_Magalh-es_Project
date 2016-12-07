%Data Analysis for two strain toggle - following ilastik and watershedding

for t=1:12 %loop over cell trapping regions
    
    trap = sprintf('%02d',t); %format the trap number into two digits to match image labels

    name = '6-20-16-toggle1-hallway001xy'; %file name
    total = 243; %total number of images per channel

    %select only trap to crop images
    img = imread(strcat(name,trap,'c1t001.tif')); %import the first (t001) phase (c1) image for current trap (xy)
    figure(1);
    imshow(img, [0 4095]);
    title('Select Trap') %select the trapping region to crop all image analysis
    rect = getrect;
    close

    %crop dimensions
    x1 = round(rect(1,1));
    x2 = round(rect(1,1) + rect(1,3));
    y1 = round(rect(1,2));
    y2 = round(rect(1,2) + rect(1,4));

    for i=1:total %loop over total number of images

        data(i,1) = 6*(i-1); %imaged ever 6 minutes, first image is t=0 min
        num = sprintf('%03d',i); %format the number into three digits to match image labels

        t = imread(strcat('watershed',trap,'t',num,'.png')); %read in the watershed mask
        phase = uint16(t);
        cropph = phase(y1:y2, x1:x2); %crop mask

        yfp = imread(strcat(name,trap,'c2t',num,'.tif')); %read in the yellow image
        cropy = yfp(y1:y2, x1:x2); %crop image
        truey = cropy.*cropph; %multiply by watershed mask

        cfp = imread(strcat(name,trap,'c3t',num,'.tif')); %read in the cyan image
        cropc = cfp(y1:y2, x1:x2); %crop image
        truec = cropc.*cropph; %multiply by watershed mask

        ratio = truey./truec; %take ratio of every yellow:cyan value

        yellow = ratio>0.5; %set ratio cut off for a yellow cell - adjustable
        cyan = ratio<=0.5 & cropph==1; %set ratio cut off for a cyan cell - adjustable without counting the entire background
        y16 = uint16(yellow); %convert to 16bit 
        c16 = uint16(cyan); %convert to 16bit 

        %save all the masks to compare to the acutal images
        imwrite(cropph, strcat('phase',trap,'t',num,'v05.png'));
        imwrite(y16, strcat('yellow',trap,'t',num,'v05.png'));
        imwrite(c16, strcat('cyan',trap,'t',num,'v05.png'));
        
        %save all the data to graph
        data(i,1) = 6*(i-1); %imaged ever 6 minutes, first image is t=0 min
        data(i,2) = sum(sum(yellow)); %numnber yellow pixels
        data(i,3) = sum(sum(cyan)); %number cyan pixels
        data(i,4) = mean2(truey); %overall yellow fluorescence
        data(i,5) = mean2(truec); %overall cyan fluorescence
        data(i,6) =  sum(sum(yellow))+sum(sum(cyan)); %total number pixels - for calculating ratio

        disp(i); %to keep track of progress

    end
    
    %save data for backup
    file = strcat(name,trap,'data,mat'); 
    save(file, 'data');

    yratio(:) = data(:,2)./data(:,6); %calculate ratio of yellow cells
    cratio(:) = data(:,3)./data(:,6); %calculate ratio of cyan cells

    %graph ratio ratio v time
    figure(1);
    plot(data(:,1),yratio(:), 'g', data(:,1),cratio(:), 'b', 'LineWidth',2);hold on
    xlabel('Time(min)', 'FontSize', 12);
    ylabel('Strain Ratio', 'FontSize', 12);
    legend('Yellow Strain', 'Cyan Strain', 'FontSize',12);
    title(strcat('Trap ', trap));
    hold off
    saveas(gcf, strcat('YellowRatioGraph',trap,'.png'));



    %graph fluor over time
    figure(2);
    plot(data(:,1),data(:,4), 'g', data(:,1),data(:,5), 'b', 'LineWidth',2);
    hold on
    xlabel('Time(min)', 'FontSize', 12);
    ylabel('Overall Fluorescence (AU)', 'FontSize', 12);
    legend('YFP', 'CFP', 'FontSize',12);
    title(strcat('Trap ', trap));
    hold off
    saveas(gcf, strcat('FluorGraph',trap,'.png'));

    %graph Y fluor v ratio
    figure(3);
    scatter(yratio(:),data(:,4), 10, data(:,1), 'filled');
    hold on
    colorbar('eastoutside');
    xlabel('Ratio Yellow Cells', 'FontSize', 12);
    ylabel('YFP Fluorescence (AU)', 'FontSize', 12);
    title(strcat('Trap ', trap));
    hold off
    saveas(gcf, strcat('YellowRvFGraph',trap,'.png'));

    %graph C fluor v ratio
    figure(4);
    scatter(cratio(:),data(:,5), 10, data(:,1), 'filled');
    hold on
    colorbar('eastoutside');
    xlabel('Ratio Cyan Cells', 'FontSize', 12);
    ylabel('CFP Fluorescence (AU)', 'FontSize', 12);
    title(strcat('Trap ', trap));
    hold off
    saveas(gcf, strcat('CyanRvFGraph',trap,'.png'));

end
