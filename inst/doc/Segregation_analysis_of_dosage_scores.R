## ----set-options, echo=FALSE, cache=FALSE---------------------------------------------------------
options(width = 100)

## -------------------------------------------------------------------------------------------------
library(fitPoly)
data(scores)

## ----echo=FALSE-----------------------------------------------------------------------------------
df <- scores[scores$marker==5,]
for (col in c(4:9,11)) df[,col] <- round(df[,col],3)
head(df)

## -------------------------------------------------------------------------------------------------
# specify parental and F1 samples:
par1 <- "P540a1"
par2 <- c("P867a1","P867a2","P867b")
F1 <- levels(scores$SampleName)[substr(levels(scores$SampleName),1,1)=="K"]
chk1 <- checkF1(scores, par1, par2, F1,
                polysomic=TRUE, disomic=TRUE, mixed=FALSE, 
                ploidy=4, outfile=NA)
#some parts of the result:
knitr::kable(chk1[1:4,1:14])
knitr::kable(chk1[1:4,c(1:2, 19:22, 26)])

## -------------------------------------------------------------------------------------------------
cordos <- correctDosages(chk1, scores, par1, par2, ploidy=4,
                        polysomic=TRUE, disomic=TRUE, mixed=FALSE)
# show the nonzero shifts:
cordos[cordos$shift!=0,]

## -------------------------------------------------------------------------------------------------
#select the markers where a shift should be tried:
cordos <- cordos[cordos$shift != 0,] 
subscores <- scores[scores$MarkerName %in% as.character(cordos$MarkerName),]
chk2 <- checkF1(subscores, par1, par2, F1,
                polysomic=TRUE, disomic=TRUE, mixed=FALSE, 
                ploidy=4, outfile=NA, shiftmarkers=cordos)

## -------------------------------------------------------------------------------------------------
chk1$shift <- 0
chk <- rbind(chk1, chk2)
chk <- chk[order(chk$MarkerName),]

## ----eval=FALSE-----------------------------------------------------------------------------------
#  data("XYdat")
#  XYgeno <- combineFiles(XYdata=XYdat, scores=scores)
#  # define qall levels of 0, 0.05, 0.10 up to 1 where we want to inspect some SNPs:
#  qall.levels <- seq(0, 1, by=0.05)
#  # select six SNPs with qall values near each of these values and draw their plots:
#  chkx <- selMarkers_qall(chk,  qall.levels, mrkperlevel=6)
#  drawXYplots(dat=XYgeno, markers=chkx,
#              out="your-path-and filename",
#              genocol=get.genocol(ploidy=4), sample.groups=list(par1, par2),
#              groups.col=c("red", "blue"), ploidy=4)

## -------------------------------------------------------------------------------------------------
#select only the markers with qall_mult > 0:
chk <- chk[!is.na(chk$qall_mult) & chk$qall_mult > 0,]
#write the dosage file, applying the shifts as listed in chk:
dosages <- writeDosagefile(chk, scores, par1, par2, F1,
                           polysomic=TRUE, disomic=TRUE, mixed=FALSE,
                           ploidy=4, scorefile=NA)
knitr::kable(dosages[10:14, 1:12])

## -------------------------------------------------------------------------------------------------
cpp <- compareProbes(chk, scores, parent1=par1, parent2=par2, F1=F1,
                     polysomic=TRUE, disomic=TRUE, mixed=FALSE,
                     ploidy=4, compfile=NA, combscorefile=NA)
knitr::kable(head(cpp$compstat))
knitr::kable(cpp$combscores[1:4, 1:12])


## -------------------------------------------------------------------------------------------------
rr <- removeRedundant(compstat=cpp$compstat, combscores=cpp$combscores,
                      compfile=NA, combscorefile=NA)
knitr::kable(rr$combscores[1:4, 1:12])


