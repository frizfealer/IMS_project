MALDI_IMS_preprocessing <- function( iFilePath, oFolderPath, conPeakSelectFlag, binningFlag ) {
#usage: 
#iFilePath: input mzML file path, e.g. "D:\\2013_Bmyc_Paeni_Early_LP_nopp.mzML"
#oFolderPath: output folder path e.g. "D:\\"
#conPeakSelectFlag: 0 or 1 flag to use consensus peak picking method.
#binningFlag: 0 or 1 flag to use binning method.
	#try to create outputFolder first
	dir.create(oFolderPath)
	inFileName = basename( iFilePath );
	#need MALDIquant and MALDIquantForeign libraries
	library("MALDIquant")
	library("MALDIquantForeign")
	library("MassSpecWavelet")
	s <- importMzMl(iFilePath)
	##check data availability
	#check if there is any empty spectrum in the list
	writeLines( 'Checking if there is empty spectrum...' );
	any(sapply(s, isEmpty))
	writeLines( 'Correcting if a data point is not in the mass spectrum class...' );
	#fixed some data, if they are read as "mass peaks"
	#e.g., sample 609 is read as "mass peak", sample 1 is read as "mass spectrum"
	ins = vector( mode = "numeric", length = length(s) );
	for (i in 1:length(s))
	{
		ins[i] = isMassSpectrum(s[[i]]);
	}
	tIdx = which( ins == 1, 1)[1];
	for (i in 1:length(s))
	{
		if( ins[i] == 0 )
		{
			idx = which(mass(s[[i]]) %in% mass(s[[tIdx]]))
			nIntVec = vector( mode = "numeric", length = length(mass(s[[tIdx]])) );
			nIntVec[idx] = intensity(s[[i]]);
			cMeta = metaData(s[[i]]);
			s[[i]] <- createMassSpectrum( mass=mass(s[[tIdx]]), intensity=nIntVec, metaData=cMeta );
		}
	}
	isMassSpectrumList(s)
	#checking if there is NA in data
	writeLines( 'Checking if there is NA values in spectra...' );
	any(!sapply(s, isRegular))
	writeLines( 'Checking if the lengths of each spectrum are equal...' );
	#some location has smaller # m/z, because MALDIQuant remove negative m/z value
	tSta = table(sapply(s, length));
	tSta
	###### peak picking #####
	pData = list()
	pDataClass = list();
	for (i in 1:length(s))
	{
		print(i)
		s1 = s[[i]]
		# whether to do normalization or not
		#s2 <- transformIntensity(s1, method="sqrt")
		# s3 <- smoothIntensity(s2, method="MovingAverage", halfWindowSize=2)
		# writeLines( 'preprocessing and peak picking...' );
		# p <- detectPeaks(s1, method="MAD", halfWindowSize=20, SNR=3)
		# pData[[i]]=p
		#signal to noise ratio
		SNR.Th = 2;
		#
		peakInfo <- peakDetectionCWT(intensity(s1), SNR.Th=SNR.Th, nearbyPeak=TRUE);
		majorPeakInfo = peakInfo$majorPeakInfo
		betterPeakInfo <- tuneInPeakInfo(intensity(s1), majorPeakInfo)
		pData[[i]] = betterPeakInfo$peakIndex
		pData[[i]] = unique( pData[[i]]);
		pDataClass[[i]] <- createMassPeaks( mass=mass(s1)[ pData[[i]] ], intensity=intensity(s1)[ pData[[i]] ], 
		metaData=list(name="peaks from CWT") )
		#plotPeak(intensity(s1), betterPeakInfo$peakIndex, range=c(1,length(intensity(s1))), main=paste('identitifed peak for surfactin LN') );
	}
	##### construct intensity Matrix and indicator Matrix #####
	dataMatrix = intensityMatrix( pDataClass, s )
	dataMatrix = t(dataMatrix)
	peakMZVec = as.double(row.names(dataMatrix));
	#construct this m/z list for comparison
	peakMZVec2 = c();
	for ( i in 1:length(s) )
	{
		peakMZVec2 = c( peakMZVec2, mass(pDataClass[[i]]) )
	}
	peakMZVec2 = sort(unique(peakMZVec2));
	#indMatrix, and indicator matrix (0 or 1) represents appearance.
	indMatrix = matrix(nrow = dim(dataMatrix)[1], ncol = dim(dataMatrix)[2]); 
	#dataMatrix2 = matrix(nrow = length(peakMZVec2), ncol = dim(dataMatrix)[2]); 
	for ( i in 1:length(s) )
	{
		idx = which( peakMZVec2 %in% mass(pDataClass[[i]]) )
		indMatrix[,i]=0;
		indMatrix[idx,i]=1;
		# dataMatrix2[,i] = 0;
		# dataMatrix2[idx,i] = intensity(pDataClass[[i]])
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
		mDataMatrix = matrix( 0, length(mPMZ), length(s) )
		mNonSDataMatrix = matrix( 0, length(mPMZ), length(s) )
		tmp = unique(indVec);
		for( i in 1:length(s) )
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
	##### get position information #####
	posMat = matrix( 0, length(s), 2 );
	for (i in 1:length(s))
	{
		target = metaData(s[[i]])$id
		target = sub("(.+)_(R[0-9]+X[0-9]+Y[0-9]+)(_.+)","\\2", target )
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