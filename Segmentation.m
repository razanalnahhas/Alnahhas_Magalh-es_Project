%% Segmentation

for ii =2:6
    trap = sprintf('%02d',ii);
    
    name = '20160929_atc8x_arab0035_crwopf_qs1_customxy';
    total = 100; % total number of images per channel

    for i=1:total

        num = sprintf('%03d',i);

        t = readIlastikFile(strcat(name,'xy',trap,'c1t',num,'_Simple Segmentation.h5')); 
        phase = uint16(t);

        imwrite(phase, strcat('phase',trap,'t',num,'.png'));

        % cleaning of mask
        m = bwareaopen(phase,300); % removes small elements
        maskdilation = imdilate(m,strel('disk',1));
        maskerosion = imerode(maskdilation,strel('disk',1));
        mask = maskerosion;

        % segmentation with watershed transform
        D = -bwdist(~mask); % using cleaner mask
        D(~mask) = -Inf;
        mm = imextendedmin(D,2);
        D2 = imimposemin(D,mm);
        Ld2 = watershed(D2);
        BW = mask;
        BW(Ld2 == 0) = 0;
        imwrite(BW, strcat('watershed',trap,'t',num,'.png'));

        disp(i);
    end
end
