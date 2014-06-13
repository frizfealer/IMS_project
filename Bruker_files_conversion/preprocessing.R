#need MALDIquant and MALDIquantForeign libraries
library("MALDIquant")
library("MALDIquantForeign")
# args <- commandArgs(trailingOnly = TRUE)
# iPath <- args[1]
iPath="D:\\Users\\YeuChern\\Dropbox\\unc\\CS\\RA\\Project_dictionaryLearning_IMS\\data\\2013 Bmyc Paeni Early LP_nopp.mzML"
s <- importMzMl(iPath)
##check data availability
#checking if there is NA in data
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

# for (i in 1:length(s))
# {
	# if( !all( mass(s[[i]]) %in% mass(s[[2]])) )
	# {
		# print(i)
	# }
# }
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

##check the preprocessing ability
## preprocessing
# s1 = s[[ti]]
# ## sqrt transform (for variance stabilization)
# s2 <- transformIntensity(s1, method="sqrt")
# ## 21 point Savitzky-Golay-Filter for smoothing spectra
# s3 <- smoothIntensity(s2, method="SavitzkyGolay", halfWindowSize=20)
# s3 <- smoothIntensity(s2, method="MovingAverage", halfWindowSize=2)
# s4 <- removeBaseline(s3, method="SNIP", iterations=100)
# p <- detectPeaks(s4, method="MAD", halfWindowSize=20, SNR=2)

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
	s2 <- transformIntensity(s1, method="sqrt")
	s3 <- smoothIntensity(s2, method="MovingAverage", halfWindowSize=2)
	s4 <- removeBaseline(s3, method="SNIP", iterations=100)
	p <- detectPeaks(s4, method="MAD", halfWindowSize=20, SNR=2.5)
	pData[[i]]=p
}
#check if the peak picked is all in the original data (No interpolation)
for (i in 1:length(s))
{
	target = mass(pData[[i]])
	if (!all( target %in% mass(s[[i]])))
	{
		print(i)
	}
}

writeLines( 'computing consensus peak at 1%...' );
cVec = vector( mode = "numeric", length = length(mass(s[[2]])) )
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

lVec = vector( mode = "integer", length = length(s) );
for (i in 1:length(s))
{
	lVec[i] = length(mass(s[[2]]));
}
lLoc = which.max(lVec);
#some location has smaller # m/z, because MALDIQuant remove negative m/z value
#location lLoc has the maximun # m/z
pMZ = mass(s[[lLoc]])[pIdx]
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
dVec = diff(pMZ)
twoDVec = dVec[1:length(dVec)-1]+dVec[2:length(dVec)]
startPos = which(twoDVec==max(twoDVec))+1 #11794
gNum = 1;
indVec = vector( mode="integer", length=length(pMZ))
indVec[startPos]=gNum
gNum = gNum + 1
tol = 0.5
currentPos = startPos+1;
ins = which( pMZ > pMZ[currentPos]+tol )
while( length(ins) != 0 )
{
	nextPos = min(ins)-1;
	indVec[currentPos:nextPos]=gNum;
	gNum = gNum + 1;
	currentPos = nextPos+1;
	ins = which( pMZ > pMZ[currentPos]+tol )
}
if( currentPos != length(pMZ) )
{
	indVec[currentPos:length(pMZ)] = gNum
	gNum = gNum + 1;
}
currentPos = startPos-1;
ins = which( pMZ < pMZ[currentPos]-tol )
while( length(ins) != 0 )
{
	nextPos = max(ins)+1;
	indVec[currentPos:nextPos]=gNum;
	gNum = gNum + 1;
	currentPos = nextPos-1;
	ins = which( pMZ < pMZ[currentPos]-tol )
}
if( currentPos != 1 )
{
	indVec[1:currentPos] = gNum
}
mPMZ = vector( mode="numeric", length=gNum )
for (i in 1:gNum)
{
	mPMZ[i]=pMZ[min( which( indVec==i ) )];
}
mPMZ = sort( mPMZ, index.return=TRUE );
mDataMatrix = matrix( 0, length(mPMZ$x), length(s) )
for( i in 1:length(s) )
{
	temp = tapply(dataMatrix[,i], indVec, max)
	mDataMatrix[,i]=temp[mPMZ$ix]
}
write.table(mDataMatrix, "dataMatrix.csv", row.names=FALSE, col.names=FALSE, sep=",")
write.table(mPMZ$x, "mzVec.csv", row.names=FALSE, col.names=FALSE, sep=",")

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
write.table(posMat, "posMat.csv", row.names=FALSE, col.names=FALSE, sep=",")
