## ----set-options, echo=FALSE, cache=FALSE---------------------------------------------------------
options(width = 100)

## -------------------------------------------------------------------------------------------------
library(fitPoly)

## -------------------------------------------------------------------------------------------------
# show part of a GenomeStudio FullDataTable file (only some relevant columns):
fname <- system.file("extdata", "FullDataTable.txt", package="fitPoly")
dat <- readDatfile(fname, check.names=FALSE)
dat[1:3, c(2, 15,16, 21,22, 27,28)]

## -------------------------------------------------------------------------------------------------
# show part of an Affymetrix AxiomGT1.summary.txt file:
fname <- system.file("extdata", "AxiomGT1.summary.txt", package="fitPoly")
dat <- readDatfile(fname, comment.char="#", check.names=FALSE)
dat[1:6, 1:4]

## -------------------------------------------------------------------------------------------------
# import a FullDataTable.txt file from GenomeStudio:
fname <- system.file("extdata", "FullDataTable.txt", package="fitPoly")
datGS <- readFullDataTable(filename=fname, out="")
head(datGS)
unique(datGS$MarkerName)
unique(datGS$SampleName)

# import an AxiomGT1.summary.txt file from Affymetrix Power Tools:
fname <- system.file("extdata", "AxiomGT1.summary.txt", package="fitPoly")
datAX <- readAxiomSummary(AXdata=fname, out="")
head(datAX)
unique(datAX$MarkerName)
unique(datAX$SampleName)

## -------------------------------------------------------------------------------------------------
fname <- system.file("extdata", "CsvAnnotationFile.v1.txt", package="fitPoly")
mrktable <- read.csv(fname, comment.char="#")
levels(datAX$MarkerName) <- 
  mrktable$customer_id[match(levels(datAX$MarkerName), mrktable$Probe.Set.ID)]

## -------------------------------------------------------------------------------------------------
fname <- system.file("extdata", "AX_sampletable.csv", package="fitPoly")
smptable <- read.csv(fname)
head(smptable)
# no split, return the original data.frame with substituted SampleNames:
datAX_unsplit <- splitNrenameSamples(dat=datAX, sampletable=smptable,
                                     SampleID="BestArray", CustomerID="SampleName",
                                     out=NA)
head(datAX_unsplit) 
# split samples based on ploidy, return a list of data.frames:
datAX_split <- splitNrenameSamples(dat=datAX, sampletable=smptable,
                                   SampleID="BestArray", CustomerID="SampleName",
                                   Ploidy="Ploidy", out=NA)
head(datAX_split[[1]])
head(datAX_split[[2]])

## -------------------------------------------------------------------------------------------------
data(XYdat)
head(XYdat)

## -------------------------------------------------------------------------------------------------
sampRmean <- tapply(XYdat$R, XYdat$SampleName, mean)
hist(sampRmean, breaks = 20)

## -------------------------------------------------------------------------------------------------
Rstats <- calcRstats(XYdat, out=NA)
quantile(Rstats$q95)

## ----eval=FALSE-----------------------------------------------------------------------------------
#  Rlevels <- seq(200, 6000, by=400)
#  sel.Rstats <- selMarkers_byR(Rstats, Rlevels, mrkperlevel=3)
#  drawXYplots(dat=XYdat, markers=sel.Rstats,
#              out="your-path-and filename", drawRthresholds=TRUE)

## -------------------------------------------------------------------------------------------------
keep <- Rstats$MarkerName[Rstats$q95 >= 1400]
dat <- XYdat[XYdat$MarkerName %in% keep,]

## -------------------------------------------------------------------------------------------------
dat.sel <- makeFitPolyFiles(XYdat, out=NA, Rquantile=0.95, marker.threshold=1400,
                            Rthreshold.param=c(0.95, 0.5, 0))
#dat.sel is in this case a list with only one useful element (the other element
#is NA); simplify:
dat.sel <- dat.sel[[1]]
head(dat.sel)

## ----eval=FALSE-----------------------------------------------------------------------------------
#  fitMarkers(ploidy=4, markers=1:3, data=dat.sel, filePrefix="A", rdaFiles=TRUE)

