function [ gW, gW0, usedElement ] = synthesizeW_fusion_terms_updated( MLEN, HEIGHT, WIDTH, type, supPerc, sparsePerc, scale, verbose)
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
        gW(i, WEnt) = (VOL_VAR*abs( randn( length(WEnt), 1 ) ) )* scale;
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
            if length(scale) >= length(usedElement)
                volInt = scale(i);
                tScale = scale(i);
            else
                volInt = scale(1);
                tScale = scale(1);
            end
            tmp = randi( 3, HEIGHT, WIDTH );
            tmp(tmp==2)=1.5; tmp(tmp==3)=2;
            target = volInt+ randn( length(find(com(j).areaIdx==1)), 1)*tScale.*tmp( com(j).areaIdx );
            target(target < 0 ) = 0;
%             tmp = randperm( length(com(j).areaIdx ), round(length(com(j).areaIdx)*0.8) );
            gW(i, com(j).areaIdx ) = target;% + tScale*tmp( com(j).areaIdx );
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
        cW = gW(i, :);
        if sparsePerc < 1
            gW(i,: ) = 0;
            for j = 1:COMMUNITY_NUM
                cIdx = find( com(j).areaIdx == 1);
                areaLen = length(cIdx)*0.8;
                areaLen = round( sqrt(areaLen) );
                [y, x] = ind2sub( [HEIGHT WIDTH],  cIdx );
                yStart = max( min(y) + randi( max(y) - min(y), 1 ) - round( areaLen / 1.5 ), min(y) );
                yEnd = min( yStart + areaLen, max(y) );
                xStart = max( min(x) + randi( max(x)-min(x), 1 ) - round( areaLen / 1.5 ), min(x) );
                xEnd = min( xStart + areaLen, max(x) );
                tmp = zeros( HEIGHT, WIDTH );
                tmp( yStart:yEnd, xStart:xEnd ) = 1;
                tmp = find( tmp == 1 );
                gW(i, tmp) = cW(1, tmp);
            end
        end
        if verbose == 1
            if cnt <= 5
                figure; imagesc( reshape( gW(i,:), HEIGHT, WIDTH ) ); colorbar;
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

