function [ gW, gW0, usedElement ] = synthesizeW( MLEN, HEIGHT, WIDTH, type, supPerc, sparsePerc, scale, verbose)
%synthesize_W synthesize W tensor for my experiment
%for different type of W, scale
%type: 'random', 'diffusion', 'thresholding'
%supNum: # dictionary element used in across all grids
%sparsePerc: the percent of element used in each grid location
%scale: the scale of weight 
if isempty( scale )
    scale = 1;
end
gW = zeros( MLEN, HEIGHT, WIDTH );
gW0 = zeros( HEIGHT, WIDTH );

cnt = 1;
if strcmp( type, 'random' ) == 1
    VOL_VAR = 1;
    usedElement = sort(randperm( MLEN, MLEN*supPerc ));
    for i = usedElement
        WEnt = sort( randperm( HEIGHT*WIDTH, HEIGHT*WIDTH*sparsePerc ) );
        gW(i, WEnt) = VOL_VAR*abs( randn( length(WEnt), 1 ) ) + scale;
        if verbose == 1
            if cnt <=5
                figure; imagesc( reshape( gW(i,:), HEIGHT, WIDTH ) ); colorbar;
                cnt = cnt + 1;
            end
        end
    end
elseif strcmp( type, 'diffusion' ) == 1
    COMMUNITY_NUM = 2;
    com = genCommunityPos( COMMUNITY_NUM, HEIGHT, WIDTH );
    usedElement = sort(randperm( MLEN, MLEN*supPerc ));
    %volInt: base intensity, with rand 20 variation
    for i = usedElement
        if verbose == 1
            fprintf('%d\n', i);
        end
        for j = 1:COMMUNITY_NUM
            volInt = 4 * scale;
            target = volInt + 3*randn(1, 1);
            if target < 0
                target = 0;
            end
            tmp = rand( HEIGHT, WIDTH );
            gW(i, com(j).areaIdx) = target + 3*tmp( com(j).areaIdx );
            gW(gW<0) = 0;
        end
    end
    cnt = 1;
    for i = usedElement
%         gaussianAry = [25 50];
%         gNum = randi(2,2);
%         h = fspecial('gaussian', [gaussianAry(gNum(1)) gaussianAry(gNum(2))], 1);
        h = fspecial( 'gaussian' );
        I = imfilter( reshape( gW(i,:), HEIGHT, WIDTH ), h );
        blurNum = 9;
        for j = 1:blurNum
            I = imfilter( reshape( I, HEIGHT, WIDTH ), h );
        end
        gW(i,:) = I(:);
        tmp = com(1).areaIdx;
        for j = 2:COMMUNITY_NUM
            tmp = tmp + com(j).areaIdx;
        end
        gW(i,~tmp) = 0;
        if verbose == 1
            if cnt <= 5
                figure; imagesc( reshape( gW(i,:), HEIGHT, WIDTH ) ); colorbar;
                cnt = cnt + 1;
            end
        end
    end
elseif strcmp( type, 'thresholding' ) == 1
    COMMUNITY_NUM = 2;
    [com] = genCommunityPos( COMMUNITY_NUM, HEIGHT, WIDTH );
    usedElement = sort(randperm( MLEN, MLEN*supPerc ));
    %volInt: base intensity, with rand 20 variation
    for i = usedElement
        if verbose == 1
            fprintf('%d\n', i);
        end
        for j = 1:COMMUNITY_NUM
            volInt = 4 * scale;
            target = volInt + 3*randn(1, 1);
            tmp = rand( HEIGHT, WIDTH );
            gW(i, com(j).areaIdx) = target + 3*tmp( com(j).areaIdx );
            gW(gW<0) = 0;
        end
    end
    cnt = 1;
    for i = usedElement
%         gaussianAry = [25 50];
%         gNum = randi(2,2);
        h = fspecial( 'gaussian' );
        I= imfilter( reshape( gW(i, :), HEIGHT, WIDTH ), h );
        blurNum = 9;
        for j = 1:blurNum
            I= imfilter( reshape( I, HEIGHT, WIDTH ), h );
        end
        gW(i,:) = I(:);
        tmp = com(1).areaIdx;
        for j = 2:COMMUNITY_NUM
            tmp = tmp + com(j).areaIdx;
        end
        gW(i,~tmp) = 0;
        QUAN_LEVEL = 15;
        [~, center] = hist( gW(i,:), QUAN_LEVEL );
        for z = 1:HEIGHT*WIDTH
            [~, minC] = min( abs( gW(i, z) - center ) );
            gW(i, z) = center(minC);
        end
        gW(i, ~tmp) = 0;
        if verbose == 1
            if cnt <= 5
                figure; imagesc( reshape( gW(i,:), HEIGHT, WIDTH ) ); colorbar;
                cnt = cnt + 1;
            end
        end
    end
    
end


end

