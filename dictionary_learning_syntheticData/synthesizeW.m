function [ gW, gW0 ] = synthesizeW( SLEN, MLEN, HEIGHT, WIDTH, type, scale, verbose)
%synthesize_W synthesize W tensor for my experiment
%for different type of W, scale
if isempty( scale )
    scale = 1;
end
if strcmp( type, 'random' ) == 1
    spsVal = 1.5;
    volInt = 2 * scale;
    gW = volInt*abs( randn( MLEN, HEIGHT, WIDTH ) );
%     gW = gW .* ( gW >= spsVal*volInt);
elseif strcmp( type, 'diffusion' ) == 1
    COMMUNITY_NUM = 2;
    com = genCommunityPos( COMMUNITY_NUM, HEIGHT, WIDTH );
    gW = zeros( MLEN, HEIGHT, WIDTH );
    %volInt: base intensity, with rand 20 variation
    for i = 1:MLEN
        if verbose == 1
            fprintf('%d\n', i);
        end
        for j = 1:COMMUNITY_NUM
            volInt = 4.2 * scale;
            target = volInt + 2*randn(1, 1);
            tmp = rand( HEIGHT, WIDTH );
            gW(i, com(j).areaIdx) = target + 3*tmp( com(j).areaIdx );
            gW(gW<0) = 0;
        end
    end
    for i = 1:MLEN
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
            if i <= 5
                figure; imagesc( reshape( gW(i,:), HEIGHT, WIDTH ) ); colorbar;
            end
        end
    end
elseif strcmp( type, 'thresholding' ) == 1
    COMMUNITY_NUM = 2;
    [com] = genCommunityPos( COMMUNITY_NUM, HEIGHT, WIDTH );
    gW = zeros( MLEN, HEIGHT, WIDTH );
    %volInt: base intensity, with rand 20 variation
    for i = 1:MLEN
        if verbose == 1
            fprintf('%d\n', i);
        end
        for j = 1:COMMUNITY_NUM
            volInt = 4.2 * scale;
            target = volInt + 2*randn(1, 1);
            tmp = rand( HEIGHT, WIDTH );
            gW(i, com(j).areaIdx) = target + 3*tmp( com(j).areaIdx );
            gW(gW<0) = 0;
        end
    end
    for i = 1:MLEN
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
            if i <= 5
                figure; imagesc( reshape( gW(i,:), HEIGHT, WIDTH ) ); colorbar;
            end
        end
    end
    
end
gW0 = ones( HEIGHT, WIDTH ) * 0;


end

