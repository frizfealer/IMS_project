MALDI_IMS_preprocessing_lessMemory <- function( iFilePath, iFolderPath, oFolderPath, conPeakSelectFlag, binningFlag ) {
	library("MALDIquant")
	library("MALDIquantForeign")
	library("mzR");
	#iFolderPath = "D:\\data_Beth\\2011\\2011-07-24 ES370(LS10) + ES3 48 h at 30 deg on 0-1X, spots_RP\\2011-07-24 ES370 + ES3_RP";
	#inFPath = "D:\\data_Beth\\2011\\mzML\\2011-07-24 ES370\ + ES3_RP.mzML";
	#e.g. oFolderPath = "D:\\"
	#testing
	# inFPath = "D:\\2013 Bmyc Paeni Early LP_nopp.mzML";
	dir.create(oFolderPath);
	inFileName = basename( iFilePath );
	nameList = list.dirs(iFolderPath, full.names = FALSE, recursive = FALSE);
	inF = openMSfile(iFilePath);
	samNum = runInfo(inF)$scanCount;

	#check if the # peaks in each sample are the same
	ins = vector( mode="integer", length = samNum );
	for (i in 1:samNum)
	{
		# print(i)
		ins[i] = peaksCount( inF, i );
	}
	table(ins);
	#get the raw data of each sample and immediately converting it to peak list

	pData = list()
	pDataClass = list();
	for (i in 1:samNum)
	{
		print(i)
		cSpec <- peaks( inF, i);
		cMLen = dim( cSpec )[1];
		cMeta = list();
		cMeta$id=nameList[i];
		#remove line with mass <= 1
		sIdx = min( which( cSpec[1:cMLen, 1]>1 ) );
		#retain lines with intensities that are positive
		tIdx = setdiff( sIdx:cMLen, which( cSpec[1:cMLen, 2] < 0 ) ) 
		s1 <- createMassSpectrum( mass=cSpec[tIdx, 1], intensity=cSpec[tIdx,2], metaData=cMeta );
		SNR.Th = 2;
		peakInfo <- peakDetectionCWT(intensity(s1), SNR.Th=SNR.Th, nearbyPeak=TRUE);
		majorPeakInfo = peakInfo$majorPeakInfo
		betterPeakInfo <- tuneInPeakInfo(intensity(s1), majorPeakInfo)
		pData[[i]] = betterPeakInfo$peakIndex
		pData[[i]] = unique( pData[[i]]);
		pDataClass[[i]] <- createMassPeaks( mass=mass(s1)[ pData[[i]] ], intensity=intensity(s1)[ pData[[i]] ], 
		metaData=list(name="peaks from CWT") )
	}
	##### construct intensity Matrix and indicator Matrix #####
	#construct this m/z list for comparison
	peakMZVec = c();
	for ( i in 1:samNum )
	{
		peakMZVec = c( peakMZVec, mass(pDataClass[[i]]) )
	}
	peakMZVec = sort(unique(peakMZVec));
	#indMatrix, and indicator matrix (0 or 1) represents appearance.
	indMatrix = matrix(nrow = length(peakMZVec), ncol = samNum); 
	dataMatrix = matrix(nrow = length(peakMZVec), ncol = samNum); 
	for ( i in 1:samNum )
	{
		idx = which( peakMZVec %in% mass(pDataClass[[i]]) )
		uIdx = which( !( 1:length(peakMZVec) %in% idx ) );
		indMatrix[,i]=0;
		indMatrix[idx,i]=1;
		dataMatrix[,i] = 0;
		dataMatrix[idx,i] = intensity(pDataClass[[i]]);
		cSpec <- peaks( inF, i);
		tIdx = which( cSpec[,1] %in% peakMZVec[uIdx] );
		dataMatrix[uIdx,i] = cSpec[tIdx,2];
	}
	##### binning, outputs are mDataMatrix and mPMZ ######
	if( binningFlag == 1 )
	{
		writeLines( 'binning...' );
		ins = rowSums( dataMatrix );
		dOrd = order( ins, decreasing=TRUE );
		tol = 5e-4;
		dPMZ = peakMZVec[dOrd];
		mPMZ = vector( mode="integer", length = length(dPMZ) );
		indVec = vector( mode="integer", length = length(dPMZ) );
		test = vector( mode="integer", length = length(dPMZ) );
		curPeakID = 1
		while( length(dPMZ) != 0 )
		{
			print(curPeakID)
			lMZ = dPMZ[1];
			rMZ = dPMZ[which( dPMZ >= lMZ-lMZ*tol & dPMZ <= lMZ+lMZ*tol )];
			rMZIdx = which( peakMZVec %in% rMZ );
			indVec[rMZIdx] = curPeakID;
			mPMZ[curPeakID] = mean(rMZ);
			test[curPeakID] = length(rMZIdx);
			curPeakID = curPeakID + 1;
			rDOrdSet = which( dPMZ %in% rMZ );
			dPMZ = dPMZ[-rDOrdSet];
		}
		peakNum = max(indVec);
		tmp = mPMZ;
		mPMZ = vector( mode="integer", length = peakNum );
		mPMZ = tmp[1:peakNum];
		mPMZ = sort( mPMZ );
		mDataMatrix = matrix( 0, length(mPMZ), samNum)
		mNonSDataMatrix = matrix( 0, length(mPMZ), samNum )
		tmp = unique(indVec);
		for( i in 1:samNum )
		{
			temp = tapply(dataMatrix[,i], indVec, max)
			mDataMatrix[,i]=temp[tmp];
		}
	} else
	{
		mDataMatrix = dataMatrix;
		mPMZ = peakMZVec;
	}
	##### consensus Peak selection, outputs are mDataMatrix and mPMZ #####
	if( conPeakSelectFlag == 1 )
	{
		writeLines( 'computing consensus peak at 1%...' );
		PEAK_INT_THRES = 10;
		cVec = vector( mode = "numeric", length = length(mPMZ) )
		#check if pData[[i]] having the maximun set of peaks
		for (i in 1:length(s))
		{
			#tMassIdx = which( mass(s[[i]]) %in% mass(pData[[i]]) )
			iPeakMassIdx = which( mDataMatrix[,i] >= PEAK_INT_THRES );
			cVec[iPeakMassIdx] = cVec[iPeakMassIdx]+1;
		}
		#consensus peak strategy
		thresPrec = 0.01;
		tNum = round( length(s)*thresPrec )
		pIdx = which(cVec>=tNum)

		mPMZ = mPMZ[pIdx]
		mDataMatrix = mDataMatrix[pIdx,]
	} 
	##### output to folder path #####
	oFileName = paste( inFileName, "_data.csv", sep="" );
	write.table(mDataMatrix, paste( oFolderPath, oFileName, sep="" ), row.names=FALSE, col.names=FALSE, sep=",")
	oFileName = paste( inFileName, "_mz.csv", sep="" );
	write.table(mPMZ, paste( oFolderPath, oFileName, sep="" ), row.names=FALSE, col.names=FALSE, sep=",")
	#get position information
	#use the nameList vector extracted from the MALDI-IMS directory
	posMat = matrix( 0, samNum, 2 );
	for (i in 1:samNum)
	{
		target = nameList[i];
		target = sub("(.+)_(R[0-9]+X[0-9]+Y[0-9]+)","\\2", target )
		xVal = as.numeric( sub("R[0-9]+X([0-9]+)Y[0-9]+","\\1", target ) );
		yVal = as.numeric( sub("R[0-9]+X[0-9]+Y([0-9]+)","\\1", target ) );
		posMat[i,1]= xVal;
		posMat[i, 2]=yVal;
	}
	oFileName = paste( inFileName, "_pos.csv", sep="" );
	write.table(posMat, paste( oFolderPath, oFileName, sep="" ), row.names=FALSE, col.names=FALSE, sep=",")
	
	oFileName = paste( inFileName, "_indPeak.csv", sep="" );
	write.table(indMatrix, paste( oFolderPath, oFileName, sep="" ), row.names=FALSE, col.names=FALSE, sep=",")
}