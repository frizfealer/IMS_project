function [ ] = customizedPlot( ax, x,y, aFontSize, xtickVec, xticklVec, tFontSize, tStr, xlFontSize, xlStr, ylFontSize, ylStr )
%customizedPlot Summary of this function goes here
plot( ax, x, y );
set( ax,'XTick', xtickVec );
assert( length(xtickVec) == size(xticklVec, 1) );
set( ax, 'XTickLabel', xticklVec );
xticklabel_rotate( [], 90, [], 'Fontsize', aFontSize, 'Fontname','arial', 'fontweight','bold' );
title( tStr, 'Fontsize', tFontSize ,'Fontname','arial', 'Fontweight', 'bold' );
xlabel( xlStr, 'Fontsize', xlFontSize, 'Fontname','arial' );
ylabel( ylStr, 'Fontsize', ylFontSize, 'Fontname','arial' );
set( ax, 'FontSize', aFontSize, 'Fontweight', 'bold' );

end

