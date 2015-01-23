iPath="D:\\IMS_DATA\\Surfactin MALDI Spectra\\Surfactin_RP.mzML"
inFileName = basename( iPath );
	#need MALDIquant and MALDIquantForeign libraries
	library("MALDIquant")
	library("MALDIquantForeign")
	s <- importMzMl(iPath)
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
	writeLines( 'preprocessing and peak picking...' );
	pData = list()
	s1 = s[[i]]
	s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=10)
	s[[i]] = s3;
	p <- detectPeaks(s3, method="MAD", halfWindowSize=20, SNR=2)
	pData[[i]]=p
	
	outputPath = "D:\\"
	oFileName = paste( inFileName, "_data.csv", sep="" );
	write.table( intensity(p), paste( outputPath, oFileName, sep="" ), row.names=FALSE, col.names=FALSE, sep=",")
	oFileName = paste( inFileName, "_mz.csv", sep="" );
	write.table( mass(p), paste( outputPath, oFileName, sep="" ), row.names=FALSE, col.names=FALSE, sep=",")