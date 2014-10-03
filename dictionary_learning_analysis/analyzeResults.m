function [ sigDict, sigM, sigIdx] = analyzeResults( expRec, DTemplate, BlkDS, SpeciesM, DIonName, mzAxis, requireInfoPath, outNum, outFolder, outPicFlag, relation, refColonName )
%if outPicFlag == 0, requireInfoPath can be [], just for the significant
%Dictionary elements output
% requireInfoPath: a data structure with three field: 
% 1. bioPath: for bio file path, 2. macPath: for machine file path, 3. 
% registeredPath: for registered info. file path
% thresDiffW: threshold for the difference of W between bacteria community
% outFolder: the folder for output files
outW = expRec.outW;
% outW0 = expRec.outW0;
outD = expRec.outD;
[mLen, ~] = size(outW);
% [sLen] = size( outD, 1 );
[IHEI, IWID] = size( BlkDS.indMap );


for i = 1:mLen
    ins(i) = max(outW(i,:));
end

%choose W at least larger then 1e-3
%and sort W according to their intensities
targetW = find( ins > 1e-3 ); 
[~,idx] = sort(ins(targetW),'descend');
targetW = targetW(idx);
% covW = cov(outW(targetW(1:31),:)' );
outW = full(outW);
outW(outW<1e-3)=0;
outD = full(outD);

%choose W has discriminat distribution between bacteria colonies
%REF_AREA: reference area, usually the area where the interaction of
%bacteria conolines happened.
%detecting the compound where the interation of colonies faciliates the
%secretion
REF_AREA = refColonName;
if strcmp( relation, 'facilitate' ) == 1
    ins2 = zeros(length(targetW),1);
    for j = 1:length(targetW)
        for i = 1:BlkDS.blkNum
            bW{i} = outW(targetW(j), BlkDS.B2GMap{i});
        end
        val = 0;
        for i = 1:BlkDS.blkNum
            if i ~= REF_AREA
                val = val + mean(bW{i});
            end
        end
        ins2(j) = mean(bW{REF_AREA}) - val;
    end
    [~,idxIns2] = sort(ins2,'descend');
end
if strcmp( relation, 'hinder' ) == 1
    ins2 = zeros(length(targetW),1);
    for j = 1:length(targetW)
        for i = 1:BlkDS.blkNum
            bW{i} = outW(targetW(j), BlkDS.B2GMap{i});
        end   
        ins2(j) = sum(bW{REF_AREA});
    end
    [~,idxIns2] = sort(ins2,'ascend');
end

idxIns2=idxIns2(1:outNum)';

if outPicFlag == 1
    % bioFilePath = 'D:\Users\YeuChern\GitHub\IMS_project\example\realExperiment_100831_113_333_136_26hr_0_1XLB_2\bioPic002.jpg';
    % macFilePath = 'D:\Users\YeuChern\GitHub\IMS_project\example\realExperiment_100831_113_333_136_26hr_0_1XLB_2\macPic_2.jpg';
    bioFilePath = requireInfoPath.bioPath;
    macFilePath = requireInfoPath.macPath;
    registeredPath = requireInfoPath.registeredPath;
    drawDictionaryElementProfile( outFolder, bioFilePath, macFilePath, registeredPath, targetW(idxIns2), expRec, SpeciesM, DTemplate, DIonName, mzAxis );
%     for j = idxIns2
%         
%         drawDictionaryElementProfile( outFolder, bioFilePath, macFilePath, registeredPath, targetElement, expRec, SpeciesM, DTemplate, DIonName, mzAxis );
%     %     if ins2(j) < thresDiffW
%     %         continue;
%     %     end
%         picNum = picNum + 1;
%         fprintf( 'output MW: %g, max intensity of W: %g\n', SpeciesM( targetW(j) ), full( ins( targetW(j) ) ) );
%         sigM(cntSigDict) = SpeciesM( targetW(j) ) ;
%         figure('Visible','Off'); 
%         imagesc( reshape( outW(targetW(j),:), IHEI, IWID ) ); axis off; export_fig(MSFilePath);
%         [~, ~, ~, ~ ] = image_registered( bioFilePath, macFilePath, MSFilePath, requireInfoPath.registeredPath, 0.2, 0 );
%         export_fig(mapMSFilePath);
%         delete(MSFilePath);
%         figure('Visible','Off');
%         h1 = subplot(1,2,1);
%         possIdx = find(DTemplate(:,targetW(j))~=0);
%         mzRange = (mzAxis(possIdx(1))-10):(mzAxis(possIdx(end))+10);
%         h=gca;
%         plot(h, mzRange,zeros(length(mzRange), 1));
%         ylim([0 1.1]);
%         lineHan = [];
%         textCell=[];
%         cnt = 1;
%         sigDict(:, cntSigDict) = outD(:, targetW(j) ); 
%         cntSigDict = cntSigDict + 1;
%         for i = 1:length(possIdx)
%             cColor = ColorSpecTable( mod(i, 7)+1,: );
%             cLine = LineStyleTable{ mod(i, 3)+1 };
%             if outD(possIdx(i), targetW(j) ) >=1e-2
%                 lineHan(cnt) = ...
%                     line( [mzAxis(possIdx(i)) mzAxis(possIdx(i))], [0 outD(possIdx(i), targetW(j) )], 'LineWidth',2,  'LineStyle', cLine, 'Color', cColor);
%                 textCell{cnt} = strcat('[', deblank(DIonName{targetW(j)}{i,:}) ,']: ', num2str(mzAxis(possIdx(i))) );
%                 cnt = cnt + 1;
%     %         if outD(possIdx(i), targetW(j) ) >=1e-2
%     %             annoText = strcat('[', deblank(DIonName{targetW(j)}(i,:)) ,']');
%     %             if outD(possIdx(i), targetW(j) ) >= 0.6
%     %                 dist = sqrt(length(annoText)^2/2)*0.5;
%     %                 text( mzAxis(possIdx(i))-3, outD(possIdx(i), targetW(j))-dist, annoText, 'FontSize', 9, 'Rotation',90);
%     %             else
%     %                 text( mzAxis(possIdx(i)), outD(possIdx(i), targetW(j)), annoText, 'FontSize', 9, 'Rotation', 90);
%     %             end
%     %         end
%             end
%         end
%         legend(lineHan, textCell);
%         h=title( ['Inferred M:', num2str(SpeciesM(targetW(j)))]);
%         xlabel( 'mzValue' ); ylabel( 'intensity' );
%     %     set(gca, 'ytick', sort(y(y~=0)));
%         % set(gca, 'xtick', x(label(1):label(end)));
%         h2 = subplot(1,2,2, 'align');
%         tmp = imread(mapMSFilePath); imshow(tmp);colorbar;
%         ax=get(h2,'Position'); set(h2,'position',[ax(1) ax(2)-0.05 ax(3)+0.1 ax(4)+0.1]);colorbar;
%         h = title( 'The distribution in the space.');
% 
%         xlabel( 'width' ); ylabel( 'height' );
%         axis([h1 h2],'square')
%         hhh=gcf; set(hhh,'pos', [1 41 1366 650] );
%         export_fig([ 'number_', num2str(picNum), '_M_', num2str(SpeciesM(targetW(j))), '_DictionaryElement.jpg']);
%         %hwplotprep
%     %     orient landscape;
%         %print ( '-dpdf', ['M_', num2str(round(SpeciesM(targetW(j)))), '_DictionaryElement'] );
%         delete(mapMSFilePath);
%     end
end
sigDict= zeros( length(mzAxis), outNum);
cntSigDict = 1;
sigM = zeros( outNum, 1 );
for j = idxIns2
    sigM(cntSigDict) = SpeciesM( targetW(j) ) ;
    sigDict(:, cntSigDict) = outD(:, targetW(j) );
    cntSigDict = cntSigDict + 1;
end
sigIdx = zeros( length(sigM), 1 );
for i = 1:length(sigM)
    sigIdx(i) = find( SpeciesM == sigM(i) );
end
end

