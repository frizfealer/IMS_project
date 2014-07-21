MALDI_IMS_preprocessing <- function( inputFilePath, outputPath ) {
#usage:
	
	#need MALDIquant and MALDIquantForeign libraries
	library("MALDIquant")
	library("MALDIquantForeign")
	# args <- commandArgs(trailingOnly = TRUE)
	# iPath <- args[1]
	#iPath="D:\\Users\\YeuChern\\Dropbox\\unc\\CS\\RA\\Project_dictionaryLearning_IMS\\data\\2013 Bmyc Paeni Early LP_nopp.mzML"
	s <- importMzMl(iPath)
	##check data availability
	#check if there is any empty spectrum in the list
	writeLines( 'Checking if there is empty spectrum...' );
	any(sapply(s, isEmpty))
	#checking if there is NA in data
	writeLines( 'Checking if there is strange values in spectra...' );
	any(!sapply(s, isRegular))
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
		s4 <- removeBaseline(s3, method="SNIP", iterations=100)
		# whether to use the original intensity or not
		s[[i]] = s4;
		p <- detectPeaks(s4, method="MAD", halfWindowSize=20, SNR=2)
		pData[[i]]=p
	}
	pData <- binPeaks( pData, tolerance=0.5 )

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
	pMZLen = length( pMZ );
	indVec = vector( mode="integer", length = pMZLen )
	curPeakID = 1;
	indVec[pMZLen] = curPeakID;
	curLoc = pMZLen;
	tol = 0.5;
	for( i in (pMZLen-1):1 )
	{
		if( abs( pMZ[i]-pMZ[curLoc] ) > tol )
		{
			curLoc = i;
			curPeakID = curPeakID + 1;
			indVec[i] = curPeakID;
		}else
		{
			indVec[i] = curPeakID;
		}
	} 
	mPMZ = vector( mode="numeric", length=curPeakID )
	cnt = 1;
	mPMZ[cnt] = pMZ[pMZLen];
	for ( i in (pMZLen-1):1 )
	{
		if( indVec[i]-indVec[i+1] == 1 )
		{
			cnt = cnt + 1;
			mPMZ[cnt] = pMZ[i];
		}
	} 
	mPMZ = sort( mPMZ );
	indVec = curPeakID - indVec + 1;
	mDataMatrix = matrix( 0, length(mPMZ), length(s) )
	for( i in 1:length(s) )
	{
		temp = tapply(dataMatrix[,i], indVec, max)
		mDataMatrix[,i]=temp
	}

	write.table(mDataMatrix, paste( outputPath, "dataMatrix.csv", sep="" ), row.names=FALSE, col.names=FALSE, sep=",")
	write.table(mPMZ, paste( outputPath, "mzVec.csv", sep="" ), row.names=FALSE, col.names=FALSE, sep=",")

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
	write.table(posMat, paste( outputPath, "posMat.csv", sep="" ), row.names=FALSE, col.names=FALSE, sep=",")
}