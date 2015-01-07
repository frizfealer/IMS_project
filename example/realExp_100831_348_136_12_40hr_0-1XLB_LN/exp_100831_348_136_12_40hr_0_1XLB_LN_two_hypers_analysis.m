%% analyze smallD results
inputPath = '/csbiohome01/ycharn/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat';
[ lVec, tVec, pVec, valVec ] = readResultsAndTest( '/csbiohome01/ycharn/IMS_project-master/two_hypers_grid_search_smallD_results', inputPath );
strVec = num2str([lVec;pVec]);
strVec1 = strVec(1:25, :);
strVec2 = strVec(26:50, :);
for i = 1:25
    strVec3(i,:) = ['(' test1(i,:) ', ' test2(i,:) ')'];
end
figure;customizedPlot( gca, 1:25, valVec, 12, 1:25, strVec3, 20, 'objective function values on leave-out data', 18, '(lambda, phi)', 18, 'objective function value' );
%% analyze largeD results
inputPath = '/csbiohome01/ycharn/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat';
[ lVec, tVec, pVec, valVec ] = readResultsAndTest( '/csbiohome01/ycharn/IMS_project-master/two_hypers_grid_search_largeD_results', inputPath );
strVec = num2str([lVec;pVec]);
strVec1 = strVec(1:25, :);
strVec2 = strVec(26:50, :);
for i = 1:25
    strVec3(i,:) = ['(' test1(i,:) ', ' test2(i,:) ')'];
end
figure;customizedPlot( gca, 1:25, valVec, 12, 1:25, strVec3, 20, 'objective function values on leave-out data', 18, '(lambda, phi)', 18, 'objective function value' );
