function [ ] = drawDictionaryElementProfile( outFolder, bioFilePath, macFilePath, registeredPath, targetElement, expRec, SpeciesM, DTemplate, DIonName, mzAxis )
% bioFilePath = 'D:\Users\YeuChern\GitHub\IMS_project\example\realExperiment_100831_113_333_136_26hr_0_1XLB_2\bioPic002.jpg';
% macFilePath = 'D:\Users\YeuChern\GitHub\IMS_project\example\realExperiment_100831_113_333_136_26hr_0_1XLB_2\macPic_2.jpg';
cPath = pwd;
if ~isempty( outFolder )
    cd( outFolder )
end
MSFilePath = 'MSFile.jpg';
mapMSFilePath = 'mapMSFile.jpg';
LineStyleTable={'-', '--',':' };
figure('Visible', 'Off')
ColorSpecTable=get(gca,'ColorOrder') ; ColorSpecTable(end-1, :) = [0 0 0]; ColorSpecTable(end, :)=[1 1 0]; close(gcf)
MarkerTable={'+','o','x','*', 's', 'd'};
outW = expRec.outW;
outD = expRec.outD;
mLen = size( outW, 1 );
[IHEI, IWID] = size( expRec.outW0 );
ins = zeros(mLen, 1);
for i = 1:mLen
    ins(i) = max(outW(i,:));
end
for j = 1:length(targetElement)
    fprintf( 'output element with M = %g, max intensity of element = %g\n', SpeciesM( targetElement(j) ), full( ins( targetElement(j) ) ) );
    figure('Visible','Off');
    imagesc( reshape( outW(targetElement(j),:), IHEI, IWID ) ); axis off; export_fig(MSFilePath); close(gcf);
    [h1, ~, ~, ~ ] = image_registered( bioFilePath, macFilePath, MSFilePath, registeredPath, 0.2, 0 );
    export_fig(mapMSFilePath); close(h1)
    delete(MSFilePath);
    figure('Visible','Off');
    h1 = subplot(1,2,1);
    possIdx = find(DTemplate(:,targetElement(j))==1);
    mzRange = (mzAxis(possIdx(1))-10):(mzAxis(possIdx(end))+10);
    h=gca;
    plot(h, mzRange,zeros(length(mzRange), 1));
    ylim([0 1.1]);
    lineHan = [];
    textCell=[];
    cnt = 1;
    for i = 1:length(possIdx)
        cColor = ColorSpecTable( mod(i, 7)+1,: );
        cLine = LineStyleTable{ mod(i, 3)+1 };
        if outD(possIdx(i), targetElement(j) ) >=1e-2
            lineHan(cnt) = ...
                line( [mzAxis(possIdx(i)) mzAxis(possIdx(i))], [0 outD(possIdx(i), targetElement(j) )], 'LineWidth',2,  'LineStyle', cLine, 'Color', cColor);
            textCell{cnt} = strcat('[', deblank(DIonName{targetElement(j)}{i,:}) ,']: ', num2str(mzAxis(possIdx(i))) );
            cnt = cnt + 1;
            %         if outD(possIdx(i), targetW(j) ) >=1e-2
            %             annoText = strcat('[', deblank(DIonName{targetW(j)}(i,:)) ,']');
            %             if outD(possIdx(i), targetW(j) ) >= 0.6
            %                 dist = sqrt(length(annoText)^2/2)*0.5;
            %                 text( mzAxis(possIdx(i))-3, outD(possIdx(i), targetW(j))-dist, annoText, 'FontSize', 9, 'Rotation',90);
            %             else
            %                 text( mzAxis(possIdx(i)), outD(possIdx(i), targetW(j)), annoText, 'FontSize', 9, 'Rotation', 90);
            %             end
            %         end
        end
    end
    legend(lineHan, textCell);
    h=title( ['Inferred M:', num2str(SpeciesM(targetElement(j)))]);
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
    export_fig([ 'element#_', num2str(targetElement(j)), '_M_', num2str(SpeciesM(targetElement(j))), '_DictionaryElement.jpg']); close(hhh)
    %hwplotprep
    %     orient landscape;
    %print ( '-dpdf', ['M_', num2str(round(SpeciesM(targetW(j)))), '_DictionaryElement'] );
    delete(mapMSFilePath);
end
cd(cPath);
end

