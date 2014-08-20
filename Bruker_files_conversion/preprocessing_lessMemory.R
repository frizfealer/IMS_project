MALDI_IMS_preprocessing <- function( iPath, inD, outputPath ) {
library("MALDIquant")
library("MALDIquantForeign")
library("mzR");
#inD = "D:\\data_Beth\\2011\\2011-07-24 ES370(LS10) + ES3 48 h at 30 deg on 0-1X, spots_RP\\2011-07-24 ES370 + ES3_RP";
#inFPath = "D:\\data_Beth\\2011\\mzML\\2011-07-24 ES370\ + ES3_RP.mzML";
#e.g. outputPath = "D:\\"
#testing
# inFPath = "D:\\2013 Bmyc Paeni Early LP_nopp.mzML";
inFileName = basename( iPath );
nameList = list.dirs(inD, full.names = FALSE, recursive = FALSE);
inF = openMSfile(iPath);
samNum = runInfo(inF)$scanCount;

# s <- importMzMl(inFPath);

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
for (i in 1:samNum)
{
	print(i)
	cSpec <- peaks( inF, i);
	cMLen = dim( cSpec )[1];
	cMeta = list();
	cMeta$id=nameList[i];
	#remove line with mass <= 1
	sIdx = min( which( cSpec[1:cMLen, 1]>1 ) );
	#remove lines with intensity is negative
	tIdx = setdiff( sIdx:cMLen, which( cSpec[1:cMLen, 2] < 0 ) ) 
	s1 <- createMassSpectrum( mass=cSpec[tIdx, 1], intensity=cSpec[tIdx,2], metaData=cMeta );
	s3 <- smoothIntensity(s1, method="SavitzkyGolay", halfWindowSize=10)
	s4 <- removeBaseline(s3, method="SNIP", iterations=100)
	# whether to use the original intensity or not
	p <- detectPeaks(s4, method="MAD", halfWindowSize=20, SNR=2)
	pData[[i]]=p
}
#pData <- binPeaks( pData, tolerance=0.5 )
writeLines( 'computing consensus peak at 1%...' );
cVec = vector( mode = "numeric", length = cMLen)
#check if pData[[i]] having the maximun set of peaks
for (i in 1:samNum)
{
tMassIdx = which( cSpec[1:cMLen, 1] %in% mass(pData[[i]]) )
cVec[tMassIdx] = cVec[tMassIdx]+1;
}
#consensus peak strategy
thresPrec = 0.01;
tNum = round( samNum*thresPrec )
pIdx = which(cVec>=tNum)
pMZ = cSpec[1:cMLen, 1][pIdx]
dataMatrix <- matrix( 0, length(pMZ), samNum )
for (i in 1:samNum)
{
	tMassIdx = which( pMZ %in% mass(pData[[i]]) )
	tMass = pMZ[tMassIdx]
	tIntIdx = which( mass(pData[[i]]) %in% tMass )
	dataMatrix[tMassIdx,i]=intensity(pData[[i]])[tIntIdx]
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
mDataMatrix = matrix( 0, length(mPMZ), samNum )
for( i in 1:samNum )
{
temp = tapply(dataMatrix[,i], indVec, max)
mDataMatrix[,i]=temp
}
oFileName = paste( inFileName, "_data.csv", sep="" );
write.table(mDataMatrix, paste( outputPath, oFileName, sep="" ), row.names=FALSE, col.names=FALSE, sep=",")
oFileName = paste( inFileName, "_mz.csv", sep="" );
write.table(mPMZ, paste( outputPath, oFileName, sep="" ), row.names=FALSE, col.names=FALSE, sep=",")
#get position information
posMat = matrix( 0, samNum, 2 );
for (i in 1:samNum)
{
target = metaData(pData[[i]])$id
target = sub("(.+)_(R[0-9]+X[0-9]+Y[0-9]+)","\\2", target )
xVal = as.numeric( sub("R[0-9]+X([0-9]+)Y[0-9]+","\\1", target ) );
yVal = as.numeric( sub("R[0-9]+X[0-9]+Y([0-9]+)","\\1", target ) );
posMat[i,1]= xVal;
posMat[i, 2]=yVal;
}
oFileName = paste( inFileName, "_pos.csv", sep="" );
write.table(posMat, paste( outputPath, oFileName, sep="" ), row.names=FALSE, col.names=FALSE, sep=",")
}