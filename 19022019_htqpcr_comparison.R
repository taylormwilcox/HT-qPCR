###############################################################################
###
###
###
###
###


### Assessment of data concordance between OpenArray and qPCR in the Rattlesnake 
setwd("~/HTqPCR_manuscript")
data <- read.csv("htqpcr_concensus_19022019.csv", na.strings = 'na')

require("irr") # package for Cohen's kapa 


### Bull (mitochondrial)
BULL_mt <- kappa2(data[,colnames(data) == "BULL_qpcr_binary" | colnames(data) == "BULL1_oa_binary"]) # bull
BULL_mt$raw <- agree(data[,colnames(data) == "BULL_qpcr_binary" | colnames(data) == "BULL1_oa_binary"])$value/100
BULL_nu <- kappa2(data[,colnames(data) == "BULL_qpcr_binary" | colnames(data) == "BULL2_oa_binary"])
BULL_nu$raw <- agree(data[,colnames(data) == "BULL_qpcr_binary" | colnames(data) == "BULL2_oa_binary"])$value/100
BRKT <- kappa2(data[,colnames(data) == "BRKT_qpcr_binary" | colnames(data) == "BRKT_oa_binary"]) # brook
BRKT$raw <- agree(data[,colnames(data) == "BRKT_qpcr_binary" | colnames(data) == "BRKT_oa_binary"])$value/100
RNBT <- kappa2(data[,colnames(data) == "RNBT_qpcr_binary" | colnames(data) == "RNBT_oa_binary"]) # rainbow
RNBT$raw <- agree(data[,colnames(data) == "RNBT_qpcr_binary" | colnames(data) == "RNBT_oa_binary"])$value/100
BRNT <- kappa2(data[,colnames(data) == "BRNT_qpcr_binary" | colnames(data) == "BRNT_oa_binary"]) # brown
BRNT$raw <- agree(data[,colnames(data) == "BRNT_qpcr_binary" | colnames(data) == "BRNT_oa_binary"])$value/100

concordance <- matrix(, ncol = 4, nrow = 5)
concordance[1,] <- c("BULL mt", BULL_mt$value, BULL_mt$p.value, BULL_mt$raw)
concordance[2,] <- c("BULL nu", BULL_nu$value, BULL_nu$p.value, BULL_nu$raw)
concordance[3,] <- c("BRKT", BRKT$value, BRKT$p.value, BRKT$raw)
concordance[4,] <- c("RNBT", RNBT$value, RNBT$p.value, RNBT$raw)
concordance[5,] <- c("BRNT", BRNT$value, BRNT$p.value, BRNT$raw)

concordance
write.csv(concordance, "rattlesnake_concordance.csv")

### Bull Trout data
require(dplyr)

data2 <- read.csv("bull_trout_data.csv", na.strings = 'na')
data2$BULL_mt <- ifelse(data2$BULL_mt_wells>0,1,0)
data2$BULL_nu <- ifelse(data2$BULL_nu_wells>0,1,0)
data2$BULL_qpcr <- ifelse(data2$BULL_qPCR_wells>0,1,0)


data2 <- distinct(data2, Sample, .keep_all= TRUE)   # remove some replicate analyses
nrow(data2)                                         # unique samples = 476
data2 <- subset(data2, is.na(data2$BULL_qpcr) == F) # remove some missing 
nrow(data2)                                         # samples with qPCR data too = 465

mtDNA_concord <- function(x){
  sum(x[,6] == x[,8], na.rm = T)/nrow(x)
}

nuDNA_concord <- function(x){
  sum(x[,7] == x[,8], na.rm = T)/nrow(x)
}

mtDNA_concord(subset(data2, data2$Plate == "Plate 1")) # 0.938
mtDNA_concord(subset(data2, data2$Plate == "Plate 2")) # 0.513
mtDNA_concord(subset(data2, data2$Plate == "Plate 3")) # 0.900
mtDNA_concord(subset(data2, data2$Plate == "Plate 4")) # 0.975
mtDNA_concord(subset(data2, data2$Plate == "Plate 5")) # 0.988
mtDNA_concord(subset(data2, data2$Plate == "Plate 6")) # 0.968

nuDNA_concord(subset(data2, data2$Plate == "Plate 1")) # 0.951
nuDNA_concord(subset(data2, data2$Plate == "Plate 2")) # 0.913
nuDNA_concord(subset(data2, data2$Plate == "Plate 3")) # 0.963
nuDNA_concord(subset(data2, data2$Plate == "Plate 4")) # 0.988
nuDNA_concord(subset(data2, data2$Plate == "Plate 5")) # 0.988
nuDNA_concord(subset(data2, data2$Plate == "Plate 6")) # 0.984

