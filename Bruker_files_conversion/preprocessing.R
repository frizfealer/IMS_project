MALDI_IMS_preprocessing <- function( iPath, outputPath, baseLineFlag, binningFlag ) {
#usage:
#e.g. inputFilePath = "D:\\Users\\YeuChern\\Dropbox\\unc\\CS\\RA\\Project_dictionaryLearning_IMS\\data\\2013_Bmyc_Paeni_Early_LP\\2013_Bmyc_Paeni_Early_LP_nopp.mzML"
#e.g. outputPath = "D:\\"
	inFileName = basename( iPath );
	#need MALDIquant and MALDIquantForeign libraries
	library("MALDIquant")
	library("MALDIquantForeign")
	# args <- commandArgs(trailingOnly = TRUE)
	# iPath <- args[1]
	#iPath="D:\\Users\\YeuChern\\Dropbox\\unc\\CS\\RA\\Project_dictionaryLearning_IMS\\data\\2013 Bmyc Paeni Early LP_nopp.mzML"
	#s <- importMzMl(iPath, centroided=TRUE)
	s <- importMzMl(iPath)
	##check data availability
	#check if there is any empty spectrum in the list
	writeLines( 'Checking if there is empty spectrum...' );
	any(sapply(s, isEmpty))
	writeLines( 'Correcting if all data points are in the mass spectrum class...' );
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
	writeLines( 'Checking if there is strange values in spectra...' );
	any(!sapply(s, isRegular))
	#fixed some data, if they are read as "mass peaks"
	#e.g., sample 609 is read as "mass peak", sample 1 is read as "mass spectrum"
	# p = 609;
	# idx = which(mass(s[[p]]) %in% mass(s[[1]]))
	# nIntVec = vector( mode = "numeric", length = length(mass(s[[1]])) );
	# nIntVec[idx] = intensity(s[[p]]);
	# cMeta = metaData(s[[p]]);
	# s[[p]] <- createMassSpectrum( mass=mass(s[[1]]), intensity=nIntVec, metaData=cMeta );
	# for (i in 1:length(s))
	# {
		# if ( TRUE %in% is.na(mass(s[[i]])) )
		# {
			# cat(i, ": mass na\n");
		# }
		# if ( TRUE %in% is.na(intensity(s[[i]])) )
		# {
			# cat(i, ": intensity na\n" );
		# }
	# }
	writeLines( 'Checking if the length of spectrum is equal...' );
	#some location has smaller # m/z, because MALDIQuant remove negative m/z value
	table(sapply(s, length))
	#correction
	for (i in 1:length(s))
	{
		cM = mass( s[[i]] );
		cI = intensity( s[[i]] );
		cMLen = length( cM );
		sIdx = min( which( cM>1 ) );
		s[[i]] <- createMassSpectrum( mass=cM[sIdx:cMLen], intensity=cI[sIdx:cMLen],
	metaData=metaData( s[[i]] ) )
	}
	# oDataMatrix = matrix( 0, length(mass(s[[2]])), length(s) )
	# for (i in 1:length(s))
	# {
		# tIntIdx = which( mass(s[[i]]) %in% mass(s[[2]]) );
		# oDataMatrix[tIntIdx,i]=intensity(s[[i]])[tIntIdx]
	# }

	# for (i in 1:length(s))
	# {
		# target = metaData(s[[i]])$id
		# target = sub("(.+)_(R[0-9]+X[0-9]+Y[0-9]+)(_.+)","\\2", target )
		# ##print(target)
		# if ( grepl("R00X055Y015", target) )
		# {
			# print(target)
			# ti = i;
		# }
	# }
	# ddd = read.table(file="D:\\Users\\YeuChern\\Dropbox\\unc\\CS\\RA\\Project_dictionaryLearning_IMS\\data\\2013_Bmyc_Paeni_Early_LP\\rawData_dump_from_bruker\\2013 Bmyc Paeni Early LP_0_R00X055Y015_1.txt", header = FALSE, sep = " " );

	#preprocessing for one single spectrum						
	# pData = list()
	# for (i in 1:length(s))
	# {
		# print(i)
		# s1 = s[[i]]
		# s2 <- transformIntensity(s1, method="sqrt")
		# s3 <- smoothIntensity(s2, method="MovingAverage", halfWindowSize=2)
		# spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=10)
		# s4 <- removeBaseline(s3, method="SNIP", iterations=100)
		# p <- detectPeaks(s4, method="MAD", halfWindowSize=20, SNR=2.5)
		# pData[[i]]=p
	# }
	#drawing preprocessing results
	# par(mfrow=c(2,3))
	# xlim <- range(mass(s1)) # use same xlim on all plots for better comparison
	# plot(s1, main="1: raw", sub="", xlim=xlim)
	# plot(s2, main="2: variance stabilization", sub="", xlim=xlim)
	# plot(s3, main="3: smoothing", sub="", xlim=xlim)
	# plot(s4, main="4: baseline correction", sub="", xlim=xlim)
	# plot(s4, main="5: peak detection", sub="", xlim=xlim)
	# points(p)
	# ## label top 20 peaks
	# top20 <- intensity(p) %in% sort(intensity(p), decreasing=TRUE)[1:20]
	# labelPeaks(p, index=top20, underline=TRUE)
	# plot(p, main="6: peak plot", sub="", xlim=xlim)
	# labelPeaks(p, index=top20, underline=TRUE)
	# par(mfrow=c(1,1))

	writeLines( 'preprocessing and peak picking...' );
	pData = list()
	for (i in 1:length(s))
	{
		print(i)
		s1 = s[[i]]
		# whether to do normalization or not
		#s2 <- transformIntensity(s1, method="sqrt")
		# s3 <- smoothIntensity(s2, method="MovingAverage", halfWindowSize=2)
		s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=10)
		if( baseLineFlag == 1 )
		{
			s4 <- removeBaseline(s3, method="SNIP", iterations=100)
		} 
		else
		{
			s4 <- s3
		}
		# whether to use the original intensity or not
		s[[i]] = s4;
		p <- detectPeaks(s4, method="MAD", halfWindowSize=20, SNR=2)
		pData[[i]]=p
	}
	#pData <- binPeaks( pData, tolerance=0.5 )

	writeLines( 'computing consensus peak at 1%...' );
	cVec = vector( mode = "numeric", length = length(mass(s[[1]])) )
	#check if pData[[i]] having the maximun set of peaks
	for (i in 1:length(s))
	{
		tMassIdx = which( mass(s[[i]]) %in% mass(pData[[i]]) )
		cVec[tMassIdx] = cVec[tMassIdx]+1;
	}
	#consensus peak strategy
	thresPrec = 0.01;
	tNum = round( length(s)*thresPrec )
	pIdx = which(cVec>=tNum)

	pMZ = mass(s[[1]])[pIdx]
	dataMatrix <- matrix( 0, length(pMZ), length(s) )
	for (i in 1:length(s))
	{
		tMassIdx = which( pMZ %in% mass(pData[[i]]) )
		tMass = pMZ[tMassIdx]
		tIntIdx = which( mass(s[[i]]) %in% tMass )
		dataMatrix[tMassIdx,i]=intensity(s[[i]])[tIntIdx]
	}
	#write.table(dataMatrix, "dataMatrix.csv", row.names=FALSE, col.names=FALSE, sep=",")
	#write.table(pMZ, "mzVec.csv", row.names=FALSE, col.names=FALSE, sep=",")

	#binning
	writeLines( 'binning...' );
	# pMZLen = length( pMZ );
	# indVec = vector( mode="integer", length = pMZLen )
	# curPeakID = 1;
	# indVec[pMZLen] = curPeakID;
	# curLoc = pMZLen;
	# tol = 0.5;
	# for( i in (pMZLen-1):1 )
	# {
		# if( abs( pMZ[i]-pMZ[curLoc] ) > tol )
		# {
			# curLoc = i;
			# curPeakID = curPeakID + 1;
			# indVec[i] = curPeakID;
		# }else
		# {
			# indVec[i] = curPeakID;
		# }
	# } 
	# mPMZ = vector( mode="numeric", length=curPeakID )
	# cnt = 1;
	# mPMZ[cnt] = pMZ[pMZLen];
	# for ( i in (pMZLen-1):1 )
	# {
		# if( indVec[i]-indVec[i+1] == 1 )
		# {
			# cnt = cnt + 1;
			# mPMZ[cnt] = pMZ[i];
		# }
	# } 
	# mPMZ = sort( mPMZ );
	# indVec = curPeakID - indVec + 1;
	# mDataMatrix = matrix( 0, length(mPMZ), length(s) )
	# for( i in 1:length(s) )
	# {
		# temp = tapply(dataMatrix[,i], indVec, max)
		# mDataMatrix[,i]=temp
	# }
	if( binningFlag == 1 )
	{
		ins = rowSums( dataMatrix );
		dOrd = order( ins, decreasing=TRUE );
		tol = 0.5;
		dPMZ = pMZ[dOrd];
		mPMZ = vector( mode="integer", length = length(dPMZ) );
		indVec = vector( mode="integer", length = length(dPMZ) );
		test = vector( mode="integer", length = length(dPMZ) );
		curPeakID = 1
		while( length(dPMZ) != 0 )
		{
			print(curPeakID)
			lMZ = dPMZ[1];
			rMZ = dPMZ[which( dPMZ >= lMZ-tol & dPMZ <= lMZ+tol )];
			rMZIdx = which( pMZ %in% rMZ );
			indVec[rMZIdx] = curPeakID;
			mPMZ[curPeakID] = lMZ;
			test[curPeakID] = length(rMZIdx);
			curPeakID = curPeakID + 1;
			rDOrdSet = which( dPMZ %in% rMZ );
			dPMZ = dPMZ[-rDOrdSet];
		}
		#checking, test function
		# for( i in 1:max(indVec) )
		# {
			# tmp = which( indVec == i );
			# if ( length(tmp) > 1 )
			# {
				# difference = diff(tmp) == 1;
				# if( any(difference==FALSE) )
				# {
					# print(i)
				# }
			# }
		# }
		peakNum = max(indVec);
		tmp = mPMZ;
		mPMZ = vector( mode="integer", length = peakNum );
		mPMZ = tmp[1:peakNum];
		mPMZ = sort( mPMZ );
		mDataMatrix = matrix( 0, length(mPMZ), length(s) )
		tmp = unique(indVec);
		for( i in 1:length(s) )
		{
			temp = tapply(dataMatrix[,i], indVec, max)
			mDataMatrix[,i]=temp[tmp];
		}
	} else
	{
		mDataMatrix = dataMatrix;
		mPMZ = pMZ
	}
	oFileName = paste( inFileName, "_data.csv", sep="" );
	write.table(mDataMatrix, paste( outputPath, oFileName, sep="" ), row.names=FALSE, col.names=FALSE, sep=",")
	oFileName = paste( inFileName, "_mz.csv", sep="" );
	write.table(mPMZ, paste( outputPath, oFileName, sep="" ), row.names=FALSE, col.names=FALSE, sep=",")

	#get position information
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
	write.table(posMat, paste( outputPath, oFileName, sep="" ), row.names=FALSE, col.names=FALSE, sep=",")
}