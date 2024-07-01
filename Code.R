library(here)
library(tidyverse)
library(limma)
library(edgeR)
library(Homo.sapiens)
library(MASS)
library(reshape2)
library(reshape)
library(RColorBrewer)
library(ggVennDiagram)
library(gplots)

# extracting data from files--------
sample1 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040918_BLOOD_36581.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample1) <- c("EntrezID", "Sample1")

sample2 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040919_BLOOD_36683.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample2) <- c("EntrezID", "Sample2")

sample3 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040920_BLOOD_36828.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample3) <- c("EntrezID", "Sample3")

sample4 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040921_BLOOD_47934.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample4) <- c("EntrezID", "Sample4")

sample5 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040922_BLOOD_58029.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample5) <- c("EntrezID", "Sample5")

sample6 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040923_BLOOD_68172.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample6) <- c("EntrezID", "Sample6")

sample7 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040940_TUMOR_31176.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample7) <- c("EntrezID", "Sample7")

sample8 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040941_TUMOR_36285.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample8) <- c("EntrezID", "Sample8")

sample9 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040945_TUMOR_36828.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample9) <- c("EntrezID", "Sample9")

sample10 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040946_TUMOR_37232.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample10) <- c("EntrezID", "Sample10")

sample11 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040947_TUMOR_37289.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample11) <- c("EntrezID", "Sample11")

sample12 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040948_TUMOR_47934.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample12) <- c("EntrezID", "Sample12")

sample13 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040949_TUMOR_68052.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample13) <- c("EntrezID", "Sample13")

sample14 <- read.csv(here("C:/Users/NIKSHEP/OneDrive/Documents/AIBD Group Project", "GSM3040950_TUMOR_78322.readCounts.txt"), sep = ' ', header = FALSE)
colnames(sample14) <- c("EntrezID", "Sample14")

# combining data
rownames(sample1) <- sample1[, 1]
sample1[, 1] <- NULL
data <- sample1
head(data)

data <- cbind(data, sample2)
data <- cbind(data, sample3)
data <- cbind(data, sample4)
data <- cbind(data, sample5)
data <- cbind(data, sample6)
data <- cbind(data, sample7)
data <- cbind(data, sample8)
data <- cbind(data, sample9)
data <- cbind(data, sample10)
data <- cbind(data, sample11)
data <- cbind(data, sample12)
data <- cbind(data, sample13)
data <- cbind(data, sample14)

data <- dplyr::select(data, c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Sample7", "Sample8", "Sample9", "Sample10", "Sample11", "Sample12", "Sample13", "Sample14"))
head(data)

# reading into DGE list
x <- DGEList(data)
head(x)

# assigning factors
group <- as.factor(c("blood", "blood", "blood", "blood", "blood", "blood", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor"))
x$samples$group <- group

# preprocessing data --------
lcpm <- cpm(x, log=TRUE)
summary(lcpm)

table(rowSums(x$counts==0)==14)
keep.exprs <- filterByExpr(x)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
samplenames <- colnames(x)

rownames <- rownames(lcpm)

sortmax <- order(lcpm[ , 1], decreasing = TRUE)

new <- as.data.frame(lcpm)
new <- arrange(new, desc(Sample1))

new["GeneID"] <- rownames(new)
new
new <- melt(new, id.vars="GeneID", value.name= "Counts", variable.name= "Sample")
head(new)
new[18000, ]
colnames(new) <- c("GeneID", "Sample", "Count")
nrow(new)

lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")

# Design matrix -----
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  BllodvsTumor = blood - tumor,
  levels = colnames(design))
contr.matrix

# voom plot -------
v <- voom(x, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)

summary(decideTests(efit, p.value = 0.05))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit, p.value = 0.05)
summary(dt)
rownames(lcpm[which(dt[,1] == -1), ])

which(tfit$p.value < 0.0005)
tfit

de.common <- which(dt[,1]!=0)
length(de.common)

vennDiagram(dt[,1], circle.col=c("turquoise", "salmon"))
write.fit(tfit, dt, file="results.txt")

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))

BllodvsTumor1 <- topTreat(tfit, coef=1, n=Inf)
which(BllodvsTumor1$adj.P.Val < 0.03)
v$genes$ENSEMBL
head(BllodvsTumor1)

BllodvsTumor <- BllodvsTumor1$ENSEMBL[1:100]
i <- which(v$genes$ENSEMBL %in% BllodvsTumor)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group,
          col=mycol, trace="none", density.info="none",
          margin=c(8,8), dendrogram="column")

which(v$genes$ENSEMBL %in% FirstVsSecond)
rownames[which(v$genes$ENSEMBL %in% FirstVsSecond)]

rownames(lcpm[which(v$genes$ENSEMBL %in% FirstVsSecond), ])
x$genes$SYMBOL[which(v$genes$ENSEMBL %in% FirstVsSecond)]

UpRegulated <- rownames(lcpm[which(dt[,1] == 1), ])
DownRegulated <- rownames(lcpm[which(dt[,1] == -1), ])

length(UpRegulated)

FirstVsSecond$ENSEMBL