setwd('C:/Users/INTEL/Documents/Atom/MC_Thesis/positives_negatives')

s <- read.table('cds')
names(s) <- NULL
cds_sizes <- as.matrix(s)

m <- read.table('totMut')
names(m) <- NULL
total_mut <- as.matrix(m)

human_exome_size <- 50000000
expec <- ( total_mut %*% t(cds_sizes) ) / human_exome_size
write.table(expec, "expec_matrix.txt", row.names=FALSE, col.names=FALSE)

p <- read.table('geneMut')
names(p) <- NULL
mut_pergene <- as.matrix(p)

chi_sqr_like <- (mut_pergene - expec)^2 / expec
chi <- data.frame(chi_sqr_like)

rownames(chi) <- c("KIRC","KIRP","SARC","KICH","COAD","READ","DLBC","PCPG","THYM","MESO","STAD","ESCA","GBM","LGG","LUAD","LUSC","UCEC","UCS","ACC","HNSC","BRCA","LAML","LIHC","CHOL","OV","TGCT","THCA","BLCA","CESC","UVM","PAAD","PRAD","SKCM")
colnames(chi) <- c("ACVR1B","ACVR2A","AJUBA","AKT1","APC","AR","ARHGAP35","ARID1A","ARID5B","ASXL1","ATM","ATR","ATRX","AXIN2","B4GALT3","BAP1","BRAF","BRCA1","BRCA2","CBFB","CCND1","CDH1","CDK12","CDKN1A","CDKN1B","CDKN2A","CDKN2C","CEBPA","CHEK2","CRIPAK","CTCF","CTNNB1","DNMT3A","EGFR","EGR3","EIF4A2","ELF3","EP300","EPHA3","EPHB6","EPPK1","ERBB4","ERCC2","EZH2","FBXW7","FGFR2","FGFR3","FLT3","FOXA1","FOXA2","GATA3","H3F3C","HGF","HIST1H1C","HIST1H2BD","IDH1","IDH2","KDM5C","KDM6A","KEAP1","KIT","KMT2B","KMT2C","KMT2D","KRAS","LIFR","LRRK2","MALAT1","MAP2K4","MAP3K1","MAPK8IP1","MECOM","MTOR","NAV3","NCOR1","NF1","NFE2L2","NFE2L3","NOTCH1","NPM1","NRAS","NSD1","PBRM1","PCBP1","PDGFRA","PHF6","PIK3CA","PIK3CG","PIK3R1","POLQ","PPP2R1A","PRX","PTEN","PTPN11","RAD21","RB1","RPL22","RPL5","RUNX1","SETBP1","SETD2","SF3B1","SIN3A","SMAD2","SMAD4","SMC1A","SMC3","SOX17","SOX9","SPOP","STAG2","STK11","TAF1","TBL1XR1","TBX3","TET2","TGFBR2","TLR4","TP53","TSHZ2","TSHZ3","U2AF1","USP9X","VEZF1","VHL","WT1")
write.table(chi, "chi-like_matrix.txt", col.names=NA)
