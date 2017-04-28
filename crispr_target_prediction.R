
#install libraries if not available
#source("https://bioconductor.org/biocLite.R")
#biocLite("CRISPRseek")
#biocLite("msa")
 


#library
library(CRISPRseek)
library(msa)
inputFilePath <- readDNAStringSet("test.fasta")  #readDNAStringSet
outputDir <- getwd()

results <- offTargetAnalysis(
inputFilePath, 
findgRNAs=TRUE,
gRNA.name.prefix="gRNA",
PAM.size=3, 
gRNA.size=20, 
PAM="NGG", 
gRNA.pattern= "^G(?:(?!T{3,}).)+$",  # Start with G and avoid three consecutive uracils in any position of a gRNA
#gRNA.pattern=  “[ACG]{4,}.{3}$”, #avoid uracil in the last 4 positions of the guide sequence 
min.score=0, 
topN=1000,
findgRNAsWithREcutOnly = TRUE,
REpatternFile=system.file("extdata", "NEBenzymes.fa",package="CRISPRseek"),
minREpatternSize=6, 
findPairedgRNAOnly = FALSE,
chromToSearch = "", 
outputDir = outputDir, 
overwrite = TRUE,
gRNAoutputName = "seq1gRNAs", 
format="fasta",
exportAllgRNAs="fasta",     #c("all", "fasta", "genbank", "no"),
fetchSequence=TRUE,
upstream=200, 
downstream=200,
#baseBeforegRNA=4,    
#baseAfterPAM=3,    
useScore=TRUE,
overlap.gRNA.positions=c(17, 18),    
useEfficacyFromInputSeq=TRUE,      
outputUniqueREs=TRUE,
foldgRNAs=TRUE, 
scoring.method="CFDscore",
subPAM.activity = hash( AA = 0,
                                    AC = 0,
                                    AG = 0.259259259,
                                    AT = 0,
                                    CA = 0,
                                    CC = 0,
                                    CG = 0.107142857,
                                    CT = 0,
                                    GA = 0.069444444,
                                    GC = 0.022222222,
                                    GG = 1,
                                    GT = 0.016129032,
                                    TA = 0,
                                    TC = 0,
                                    TG = 0.038961039,
                                    TT = 0))

