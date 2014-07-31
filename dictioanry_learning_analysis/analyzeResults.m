function [ ] = analyzeResults( expRec, DTemplate, BlkDS, SpeciesM, DIonName, mzAxis, requireInfoPath, thresDiffW, outFolder )
%UNTITLED3 Summary of this function goes here
% requireInfoPath: a data structure with three field: 
% 1. bioPath: for bio file path, 2. macPath: for machine file path, 3. 
% registeredPath: for registered info. file path
% thresDiffW: threshold for the difference of W between bacteria community
% outFolder: the folder for output files
cd( outFolder );
outW = expRec.outW;
% outW0 = expRec.outW0;
outD = expRec.outD;
[mLen, ~] = size(outW);
% [sLen] = size( outD, 1 );
[IHEI, IWID] = size( BlkDS.indMap );


for i = 1:mLen
    ins(i) = max(outW(i,:));
end

targetW = find( ins > 1e-3 ); 
[~,idx] = sort(ins(targetW),'descend');
targetW = targetW(idx);
covW = cov(outW(targetW(1:31),:)' );
outW = full(outW);
outD = full(outD);

ins2 = zeros(length(targetW),1);
for j = 1:length(targetW)
    for i = 1:BlkDS.blkNum
        bW{i} = outW(targetW(j), BlkDS.B2GMap{i});
    end
    tmp = [];
    for i = 1:BlkDS.blkNum-1
        for z = (i+1):BlkDS.blkNum
            tmp = [tmp abs(mean(bW{i})-mean(bW{z}))];
        end
    end
    ins2(j) = max(tmp);
end

% bioFilePath = 'D:\Users\YeuChern\GitHub\IMS_project\example\realExperiment_100831_113_333_136_26hr_0_1XLB_2\bioPic002.jpg';
% macFilePath = 'D:\Users\YeuChern\GitHub\IMS_project\example\realExperiment_100831_113_333_136_26hr_0_1XLB_2\macPic_2.jpg';
bioFilePath = requireInfoPath.bioPath;
macFilePath = requireInfoPath.macPath;
MSFilePath = 'MSFile.jpg';
mapMSFilePath = 'mapMSFile.jpg';
for j = 1:length(targetW)
    if ins2(j) < thresDiffW
        continue;
    end
    fprintf( 'output MW: %g\n', SpeciesM( targetW(j) ) );
    figure('Visible','Off'); imagesc( reshape( outW(targetW(j),:), IHEI, IWID ) ); axis off; export_fig(MSFilePath);
    [~, ~, ~, ~ ] = image_registered( bioFilePath, macFilePath, MSFilePath, requireInfoPath.registeredPath, 0.2, 0 );
    export_fig(mapMSFilePath);
    delete(MSFilePath);
    figure('Visible','Off');
    h1 = subplot(1,2,1);
    possIdx = find(DTemplate(:,targetW(j))~=0);
    mzRange = (mzAxis(possIdx(1))-10):(mzAxis(possIdx(end))+10);
    h=gca;
    plot(h, mzRange,zeros(length(mzRange), 1));
    ylim([0 1.1]);
    for i = 1:length(possIdx)
        line( [mzAxis(possIdx(i)) mzAxis(possIdx(i))], [0 outD(possIdx(i), targetW(j) )], 'LineWidth',2 );
        if outD(possIdx(i), targetW(j) ) >=1e-3
            annoText = strcat('[', deblank(DIonName{targetW(j)}(i,:)) ,']');
            if outD(possIdx(i), targetW(j) ) >= 0.6
                dist = sqrt(length(annoText)^2/2)*0.5;
                text( mzAxis(possIdx(i))-3, outD(possIdx(i), targetW(j))-dist, annoText, 'FontSize', 9, 'Rotation',90);
            else
                text( mzAxis(possIdx(i)), outD(possIdx(i), targetW(j)), annoText, 'FontSize', 9, 'Rotation', 90);
            end
        end
    end
    h=title( ['Inferred M:', num2str(SpeciesM(targetW(j)))]);
    xlabel( 'mzValue' ); ylabel( 'intensity' );
%     set(gca, 'ytick', sort(y(y~=0)));
    % set(gca, 'xtick', x(label(1):label(end)));
    h2 = subplot(1,2,2, 'align');
    tmp = imread(mapMSFilePath); imshow(tmp);colorbar;
    ax=get(h2,'Position'); set(h2,'position',[ax(1) ax(2)-0.05 ax(3)+0.1 ax(4)+0.1]);colorbar;
    h = title( 'The distribution in the space.');

    xlabel( 'width' ); ylabel( 'height' );
    axis([h1 h2],'square')
    hhh=gcf; set(hhh,'pos', [1 41 1366 650] );
    export_fig(['M_', num2str(round(SpeciesM(targetW(j)))), '_DictionaryElement.png']);
    %hwplotprep
%     orient landscape;
    %print ( '-dpdf', ['M_', num2str(round(SpeciesM(targetW(j)))), '_DictionaryElement'] );
    delete(mapMSFilePath);
end

end

