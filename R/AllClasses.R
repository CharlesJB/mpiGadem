setClass("parameters", representation( #GADEM parameters
numWordGroup="numeric",
numTop3mer="numeric",
verbose="numeric",
numTop4mer="numeric",
numTop5mer="numeric",
numGeneration="numeric",
populationSize="numeric",
pValue="numeric",
eValue="numeric",
extTrim="numeric",
minSpaceWidth="numeric",
maxSpaceWidth="numeric",
useChIPscore="numeric",
numEM="numeric",
fEM="numeric",
widthWt="numeric",
fullScan="numeric",
slideWinPWM="numeric",
stopCriterion="numeric",
numBackgSets="numeric",
weightType="numeric",
bFileName="character",
Spwm="character",
nSequences="numeric",
maskR="numeric",
nmotifs="numeric"
))
	
setClass("align", representation(
#GADEM output :
#[10bp flanking--MOTIF INSTANCE--10bp flanking] [strand] [seqID] [pos] [p-value]
#>chr1_1102_1302	tttttcatgaCTTTTCCAGGAAAGGGGTGGGGAATTCCTAGAACTGAGGGTccctcccctt	+	10	147	4.765132e-12
seq="character", #sequence like atgAACCTgga
chr="character", #sequence chr	
start="numeric", #sequence start
end="numeric", #sequence end
strand="character", #motif strand
seqID="numeric", #read number
pos="numeric", #motif position
pval="numeric", #pvalue
fastaHeader="numeric" #Fasta accession	
))

setClass("motif", representation(
# Better to call it motif
# Contains PWM, motif consensus, motif length and all aligned sequences for a specific motif
pwm="matrix",
consensus="character",
alignList="list",
name="character"
))

setClass("gadem", representation(
#general GADEM class; contains all gadem output information
parameters="list",
motifList="list"#PWM, motif, p-value,...
))

#gadem@motifs[[1]]@align[[1]]@seq




 
