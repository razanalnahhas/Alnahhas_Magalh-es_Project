%% Analysis after Matlab segmentation (Watershed transform) - single strain toggle

for ii = 1:12;
    trap = sprintf('%02d',ii);
    
    name = '20160929_atc8x_arab0035_crwopf_qs1_customxy';
    total = 100; %total number of images per channel
    
    %select only trap to crop images
    img = imread(strcat(name,trap,'c1t001.tif'));
    figure(1);
    imshow(img, [0 4095]);
    title('Select Trap')
    rect = getrect;
    close

    %crop dimensions
    x1 = round(rect(1,1));
    x2 = round(rect(1,1) + rect(1,3));
    y1 = round(rect(1,2));
    y2 = round(rect(1,2) + rect(1,4));

    for i=1:total

        data(i,1) = 6*(i-1); %imaged ever 6 minutes, first image is t=0 min
        num = sprintf('%03d',i);

        t = imread(strcat('watershed',trap,'t',num,'.png'));
        phase = uint16(t);
        cropph = phase(y1:y2, x1:x2);

        yfp = imread(strcat(name,trap,'c3t',num,'.tif')); 
        cropy = yfp(y1:y2, x1:x2); % use the x and y values already used 
        truey = cropy.*phase;

        cfp = imread(strcat(name,trap,'c4t',num,'.tif')); 
        cropc = cfp(y1:y2, x1:x2);
        truec = cropc.*phase;

        mcherry = imread(strcat(name,trap,'c2t',num,'.tif'));
        cropm = mcherry(y1:y2, x1:x2);

        ratio = truey./truec;
        s = size(ratio);
        percent(i) = (sum(sum(ratio > 0)))/(s(1,1)*s(1,2)); 

        imwrite(truey, strcat('yellow',trap,'t',num,'.png'));
        imwrite(truec, strcat('cyan',trap,'t',num,'.png'));

        data(i,4) = mean2(truey);
        data(i,5) = mean2(truec);
        data(i,7) = mean2(cropm);
        data(i,8) = mean2(truey)+mean2(truec); % total number of pixels

        file = strcat(name,'data.mat');
        save(file, 'data');

        disp(i);

    end

    % translate mcherry intensity into inducer concentration
    ma = max(data(:,7));
    m = data(:,7)./ma;
    truem = m.*maxinducer; 

    % ratios
    yyratio(:) = data(:,4)./data(:,8); % ratio YFP/total (percentage)
    ccratio(:) = data(:,5)./data(:,8); % ratio CFP/total (percentage)
    ratioyc(:) = data(:,4)./data(:,5); % ratio YFP/CFP

    % ratio v time
    figure; 
    plot(data(:,1),ratioyc(:),'r','LineWidth',2);hold on
    xlabel('Time(min)', 'FontSize', 12);
    ylabel('Ratio YFP/CFP', 'FontSize', 12);
    hold off
    saveas(gcf, strcat('RatioTimeGraph',trap,'.png'));

    % ratio v inducer 
    figure; 
    plot(truem(:),ratioyc(:),'r','LineWidth',2);hold on
    xlabel('Inducer concentration (ng/ml)', 'FontSize', 12);
    ylabel('Ratio YFP/CFP', 'FontSize', 12);
    hold off
    saveas(gcf, strcat('RatioInducer',trap,'.png'));

    %graph fluor v time
    figure;
    plot(data(:,1),data(:,4),'y',data(:,1),data(:,5),'c','LineWidth',2);
    hold on
    xlabel('Time(min)', 'FontSize', 12);
    ylabel('Overall Fluorescence (AU)', 'FontSize', 12);
    legend('YFP', 'CFP', 'FontSize',12);
    hold off
    saveas(gcf, strcat('FluorTime',trap,'.png'));

    %graph fluor v inducer 
    figure;
    plot(truem(:),data(:,4),'y',truem(:),data(:,5),'c','LineWidth',2);
    hold on
    xlabel('Inducer concentration (ng/ml)', 'FontSize', 12);
    ylabel('Fluorescence (AU)', 'FontSize', 12);
    legend('YFP','CFP','FontSize',12);
    hold off
    saveas(gcf, strcat('FluorInducer',trap,'.png'));

    % percentage of 2-color cells v time
    figure; 
    plot(data(:,1), percent(:),'g','LineWidth',2);
    hold on
    xlabel('Time (min)', 'FontSize', 12);
    ylabel('Percentage of cells with both colors', 'FontSize', 12);
    ylim([0 1]);
    hold off 
    saveas(gcf, strcat('2colorTime',trap,'.png'));

    % percentage of 2-color cells v inducer
    figure; 
    plot(truem(:), percent(:),'g','LineWidth',2);
    hold on
    xlabel('Inducer concentration (ng/ml)', 'FontSize', 12);
    ylabel('Percentage of cells with both colors', 'FontSize', 12);
    ylim([0 1]);
    hold off
    saveas(gcf, strcat('2colorInducer',trap,'.png'));
end

