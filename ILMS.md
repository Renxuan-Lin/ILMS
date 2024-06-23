# ILMS

This repository contains codes necessary for analysis of the ILMS.

## Consensus clusters

```R
rm(list = ls())
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(GSVA)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)


#####input tcga data----
load("TCGA_data/PRAD/outputdata/PRAD_lnc_mRNA_fpkm_clin_498.rda")
###immune ssGSEA analysis----
immunegenes <- read.csv(file = "other_data/immunegenes.csv",sep=",",header = T)
list<- split(as.matrix(immunegenes)[,1], immunegenes[,2])

gsva_data <- gsva(as.matrix(mrna_fpkm_T),list, method = "ssgsea")
ss <- gsva_data

##Consensus clusters----
library(ConsensusClusterPlus)
# ss1 <- sweep(ss,1,apply(ss,1,median))
# boxplot(ss1[,1:20])
dir.create('figure_final/ConsensusCluster/')
results = ConsensusClusterPlus(as.matrix(ss),
                               maxK=9,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               tmyPal = c('navy','darkred'),
                               title='figure_final/ConsensusCluster/',
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot="pdf")

icl <- calcICL(results,title = 'figure_final/ConsensusCluster/',plot = 'pdf')

## PAC = Proportion of ambiguous clustering 
Kvec = 2:9
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK

PAC <- as.data.frame(PAC)
PAC$K <- 2:9
library(ggplot2)
ggplot(PAC,aes(factor(K),PAC,group=1))+
  geom_line()+
  theme_bw(base_rect_size = 1.5)+
  geom_point(size=4,shape=21,color='#ffbc78',fill='#ff9897')+
  ggtitle('Proportion of ambiguous clustering')+
  xlab('Cluster number K')+ylab(NULL)+
  theme(axis.text = element_text(size=12),
        plot.title = element_text(hjust=0.5),
        axis.title = element_text(size=13))
ggsave(filename = 'figure_final/ConsensusCluster/PAC.pdf',width = 3.8,height = 4)


clusterNum=2      
cluster=results[[clusterNum]][["consensusClass"]]

sub <- data.frame(Sample=names(cluster),Cluster=cluster)
sub$Cluster <- paste0('C',sub$Cluster)
table(sub$Cluster)

head(sub)

my <- results[[2]][["ml"]]
library(pheatmap)
rownames(my) <- sub$Sample
colnames(my) <- sub$Sample
pheatmap(1-my,show_colnames = F,show_rownames = F,
         treeheight_row = 20,treeheight_col = 20,
         clustering_method = 'complete',
         color = colorRampPalette(c("white","#EA8379"))(50),
         annotation_names_row = F,annotation_names_col = F,
         annotation_row = data.frame(Cluster=sub$Cluster,row.names = sub$Sample),
         annotation_col = data.frame(Cluster=sub$Cluster,row.names = sub$Sample),
         annotation_colors = list(Cluster=c('C1'="#299D8F",'C2'="#D87659")))
library(export)
graph2pdf(file='figure_final/ConsensusCluster/cluster_for_2.pdf',width=5.5,height=4.5)

# -------------------------------------------------------------------------
load("output_final/0/ssgsva_cluster.Rda")
ss2 <- merge(sub,t(ss),by.x=1,by.y=0)
ss2 <- pivot_longer(ss2,3:30,names_to = 'cell',values_to = 'value')

ggplot(ss2,aes(cell,value,fill=Cluster))+
  geom_boxplot(outlier.colour = NA)+
  stat_compare_means(label = 'p.signif')+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))
library(export)
graph2pdf(file='figure_final/ConsensusCluster/boxplot_immune_for_2.pdf',width=20,height=20)


###save data for cluster----
save(sub,ss,file = "output_final/0/ssgsva_cluster.Rda")

# -------------------------------------------------------------------------
#prepare clinical data
identical(rownames(sub),rownames(clin))
TCGA_clin <- cbind(sub,clin)##合并临床信息

##取出自己想看到的临床信息
colnames(TCGA_clin)
my <- TCGA_clin[,c(2,3,4,5,7,8,9,11)]
my$age <- ifelse(my$age>65,'>65','≤65')
table(my$age)

my <- my[order(my$Cluster,my$age,
               my$gleason_score,my$OS,my$BCR),]

ee <- as.data.frame(ss)[,rownames(my)]
ee <- t(scale(t(ee)))
table(my$Cluster)

# -------------------------------------------------------------------------

my[is.na(my)] <- 'NA'
my$age <- factor(my$age,levels = c('≤65','>65'))
my$gleason_score <- factor(my$gleason_score,levels = c('6','7','8','9',"10",'NA'))
my$Cluster <- factor(my$Cluster)
my$pT <- factor(my$pT)
my$pN <- factor(my$pN)
my$cM <- factor(my$cM)
my$BCR<- factor(my$BCR)
my$OS <-  factor(my$OS)

# -------------------------------------------------------------------------
table(my$Cluster)
Cluster <- c("#299D8F","#D87659")
names(Cluster) <- levels(my$Cluster)

age <- c(pal_nejm(alpha = 0.9)(8)[3],'#CF4E27')
names(age) <- levels(my$age)

table(my$gleason_score)
gleason_score <- c('cornsilk','paleturquoise','goldenrod','firebrick','red','white')
names(gleason_score) <- levels(my$gleason_score)

table(my$BCR)
BCR <- c('#E0864A','rosybrown')
names(BCR) <- levels(my$BCR)

table(my$OS)
OS <- c('lavenderblush','slategray',"white")
names(OS) <- levels(my$OS)

table(my$pT)
pT <- c("white",pal_npg("nrc", alpha = 0.7)(6))
names(pT) <- levels(my$pT)

table(my$pN)
pN <- c(pal_igv("default",alpha = 0.7)(2),"white")
names(pN) <- levels(my$pN)

table(my$cM)
cM <- c(pal_jco("default",alpha = 0.7)(4),"white")
names(cM) <- levels(my$cM)



# -------------------------------------------------------------------------

Top = HeatmapAnnotation(Cluster=my$Cluster,
                        age=my$age,
                        gleason_score= my$gleason_score,
                        pT=my$pT,
                        pN=my$pN,
                        cM=my$cM,
                        BCR = my$BCR,
                        OS=my$OS,
                        
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,
                                                     title_gp = gpar(fontsize = 10,fontface = "bold"),
                                                     ncol=1),
                        border = T,
                        col=list(Cluster = Cluster,
                                 age = age,
                                 gleason_score =  gleason_score,
                                 pT=pT,
                                 pN=pN,
                                 cM=cM,
                                 BCR=BCR,
                                 OS=OS
                                 
                        ),
                        show_annotation_name = TRUE,
                        annotation_name_side="left",
                        annotation_name_gp = gpar(fontsize = 10))
table(sub$Cluster)
Heatmap(ee,name='Z-score',
        top_annotation = Top,
        cluster_rows = T,
        col=colorRamp2(c(-2,0,2),c("#ABD3E1" ,'white',"#FAC795")),#49b0d9
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        column_split = c(rep(1,200),rep(2,298)),
        gap = unit(1, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 10), border = T,
                                  title_gp = gpar(fontsize = 10, fontface = "bold")),
        column_gap = unit(2,'mm')
) 
library(export)
graph2pdf(file='figure_final/ConsensusCluster/cell-heatmap.pdf',width=15,height=15)


save(my,TCGA_clin,file = "output_final/0/sgsva_cluster_clin.Rda")
# -------------------------------------------------------------------------



tt <- cbind(Cluster=as.character(my$Cluster),t(ee))
tt <- as.data.frame(tt)
tt2 <- pivot_longer(tt,cols = 2:29,names_to = 'cell',values_to = 'value')
tt2$value <- as.numeric(tt2$value)

# -------------------------------------------------------------------------

source("codes/GeomSplitViolin.R")

ggplot(tt2, aes(cell,value, fill = Cluster)) + 
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75), #画4分位线
                    trim = T, #是否修剪小提琴图的密度曲线
                    linetype = "solid", #周围线的轮廓
                    color = "black", 
                    size = 0.35,
                    na.rm = T,
                    position ="identity")+ #周围线粗细
  ylab("Relative Infiltration") + xlab(NULL) +
  scale_fill_manual(values = c("#FAC795","#ABD3E1" ))+
  stat_compare_means(label = 'p.signif')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1,size=10),
        legend.position=c(0.03,0.003),legend.justification = c(0,0),
        axis.title.y = element_text(size=12))
ggsave(filename = 'figure_final/ConsensusCluster/cell-violineplot.pdf',width = 10,height=8)


estimate <- function(dat,pro){
  library(estimate)
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")   ## 注意这个platform参数
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}

pro='PRAD'###for project
dat <- mrna_fpkm_T
scores=estimate(dat,pro)

rownames(scores)=colnames(dat)
scores=scores[rownames(my),]
identical(rownames(scores),rownames(my))
data <- cbind(scores,my)
colnames(data)
class(my$Cluster)
Cluster <- c('#eeba4d','#21b6af')
a <- data
library(ggpubr)
ggplot(a,aes(x=Cluster,y=ImmuneScore,color=Cluster))+
  geom_boxplot()+
  geom_jitter(shape=16, position = position_jitter(0.2))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = c('none')
  )+
  scale_color_manual(values=c('#eeba4d','#21b6af'))+
  stat_compare_means(label = "p.signif",method="t.test",label.x = 1.5
  )+
  xlab("")
graph2pdf(file='figure_final/ConsensusCluster/estimate_for_2.pdf',width=10,height=10)
```

## immune score

```R
library(IOBR)
# MCPcounter----
im_mcpcounter <- deconvo_tme(eset = dat,
                             method = "mcpcounter"
)

# EPIC ----
im_epic <- deconvo_tme(eset = dat,
                       method = "epic",
                       arrays = F
)


# xCell ----
im_xcell <- deconvo_tme(eset = dat,
                        method = "xcell",
                        arrays = F
)


# CIBERSORT ----
im_cibersort <- deconvo_tme(eset = dat,
                            method = "cibersort",
                            arrays = F,
                            perm = 1000
)



# IPS ----
im_ips <- deconvo_tme(eset = dat,
                      method = "ips",
                      plot = F
)


# quanTIseq ----
im_quantiseq <- deconvo_tme(eset = dat,
                            method = "quantiseq",
                            scale_mrna = T
)


# ESTIMATE ----
im_estimate <- deconvo_tme(eset = dat,
                           method = "estimate"
)

# TIMER ----
dim(dat)
im_timer <- deconvo_tme(eset = dat
                        ,method = "timer"
                        ,group_list = rep("prad",dim(dat)[2])
)

save(im_cibersort,im_epic,im_ips,im_mcpcounter,
     im_quantiseq,im_timer,im_xcell,im_estimate,
     file = "output_final/1/immue/immuneinflation.Rda")


##整合数据----
library(tidyverse)
load("output_final/1/immue/immuneinflation.Rda")
load("output_final/0/ssgsva_cluster.Rda")
tme <- list(im_cibersort=im_cibersort[,-c(24:26)],
            im_epic=im_epic,
            im_ips=im_ips,
            im_mcpcounter=im_mcpcounter,
            im_quantiseq=im_quantiseq,
            im_timer=im_timer,
            im_xcell=im_xcell)
tme <- lapply(tme,function(x){
  x[,-1] <- scale(x[,-1])
  return(x)})
tme_combine <- tme$im_cibersort[,-6]%>% 
  inner_join(tme$im_epic, by="ID") %>% 
  inner_join(tme$im_ips, by="ID") %>% 
  inner_join(tme$im_mcpcounter, by= "ID") %>% 
  inner_join(tme$im_quantiseq, by= "ID") %>% 
  inner_join(tme$im_timer, by= "ID") %>% 
  inner_join(tme$im_xcell, by= "ID") %>% 
  column_to_rownames("ID")
annotation_names <- c("CIBERSORT", "EPIC", "ips", "MCPCounter",
                      "quantiSeq", "TIMER", "xCell")
annotation_cibersort <- rep(annotation_names[1], ncol(im_cibersort)-5)
annotation_epic <- rep(annotation_names[2], ncol(im_epic)-1)
annotation_ips <- rep(annotation_names[3], ncol(im_ips)-1)
annotation_MCPCounter<- rep(annotation_names[4], ncol(im_mcpcounter)-1)
annotation_quantiSeq <- rep(annotation_names[5], ncol(im_quantiseq)-1)
annotation_TIMER <- rep(annotation_names[6], ncol(im_timer)-1)
annotation_xCell <- rep(annotation_names[7], ncol(im_xcell)-1)

annotation <- c(annotation_cibersort, annotation_epic, annotation_ips,
                annotation_MCPCounter,
                annotation_quantiSeq,annotation_TIMER,
                annotation_xCell)

method=cbind(colnames(tme_combine),annotation) %>% 
  as.data.frame() %>% 
  column_to_rownames("V1") %>% as.data.frame() 

method$annotation=factor(method$annotation)
method=method$annotation
sub <- sub[order(sub$Cluster),]
table(sub$Cluster)

k <- rownames(sub)
tme_combine <- tme_combine[k,]
identical(rownames(sub),rownames(tme_combine))

tme_combine2 <- cbind(group=sub$Cluster,tme_combine)
tme_combine2$group<-factor(tme_combine2$group,levels = c("C1","C2"))
library(ComplexHeatmap)
left <- rowAnnotation(
  foo = anno_block(
    gp = gpar(fill = 2:8),
    labels =  c("CIBERSORT", "EPIC", "ips","MCPCounter",
                "quantiSeq", "TIMER", "xCell"),
    labels_gp = gpar(col = "white", fontsize = 3.5)
  ))


im <- t(tme_combine2[,-1])
group <- c("#299D8F","#D87659")

names(group) <- levels(tme_combine2$group)

Top1 = HeatmapAnnotation(group=tme_combine2$group,
                         annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,
                                                      title_gp = gpar(fontsize = 10,fontface = "bold"),
                                                      ncol=1),
                         border = T,
                         col= list(group=group),
                         show_annotation_name = TRUE,
                         annotation_name_side="left",
                         annotation_name_gp = gpar(fontsize = 10))
library(circlize)
Heatmap(im,name='Z-score',
        top_annotation = Top1,
        left_annotation = left,
        cluster_rows = F,
        col=colorRamp2(c(-2,0,2),c("#ABD3E1" ,'white',"#FAC795")),#49b0d9
        color_space = "RGB",
        cluster_columns = FALSE,
        border = T,
        row_title = NULL,
        column_labels = NULL,
        row_order=NULL,
        row_names_side = 'right',
        column_order=NULL,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 5),
        column_split = c(rep(1,200),rep(2,298)),
        row_split = list(group = method),
        gap = unit(1, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 5),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 10), border = T,
                                  title_gp = gpar(fontsize = 10, fontface = "bold")),
        column_gap = unit(2,'mm')
) 
```

## ImmuLnc

```R
### http://bio-bigdata.hrbmu.edu.cn/ImmReg/index.jsp
library(tidyverse)
b=read.table(file = "other_data/lncRNA_Pathway_sig.txt",quote = "",sep = "\t",header = T)
table(b$Cancer)
PRAD_lnc <- b[b$Cancer=="PRAD",]
# lncRES scores >0.995 and FDR <0.05
tt <-PRAD_lnc %>% 
  filter(abs(Score)>0.995, P.Adjust < 0.05)
tts <- unique(tt$lncRNA.Symbol)

##plot
library(VennDiagram)
venn.plot <- venn.diagram(
  list(lnc_COM=l1,Immulnc=tts),
  filename = "figure_final/1/COMGEN/lnccom.tiff",
  scaled=F,
  lty = 1,
  lwd = 1,
  col = "black",  ##圈的颜色
  fill = c("darkcyan", "darkorange1"),
  alpha = 0.60,
  cat.dist=0.02,
  cat.pos= -180,
  cat.col = "black",
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.07,
  cex = 0.8
)
```

## MOVIC

```
###select gene----
library(tidyverse)
load("all_gene1_data/TCGA.Rda") ###此处导入的都是log2处理后的数据
load("all_gene1_data/ILS.rda")
load("all_gene1_data/IMS.rda")
clin <- clin[,c(9,10)] %>% na.omit()
identical(rownames(clin),colnames(TCGA_16968))
colnames(clin) <- c("fustat","futime")
mrna <- TCGA_16968[IMS1,]
keep1 <- rowSums(mrna)>0
mrna <- mrna[keep1,]
lncrna <- TCGA_16968[ILS1,] 
keep2 <-rowSums(lncrna)>0
lncrna <- lncrna[keep2,]

TCGA <- list(mrna=mrna,
             lncrna=lncrna,
             clin=clin)

library(MOVICS)
tmp1 <- TCGA$mrna
elite.tmp1 <- getElites(dat     = tmp1,
                        method    = "cox",
                        surv.info = clin, # must provide survival information with 'futime' and 'fustat'
                        p.cutoff  = 0.05)
table(elite.tmp1$unicox$pvalue < 0.05)
# FALSE  TRUE 
# 42    34 
tmp2 <- TCGA$lncrna

elite.tmp2<- getElites(dat       = tmp2,
                       method    = "cox",
                       surv.info = clin, # must provide survival information with 'futime' and 'fustat'
                       p.cutoff  = 0.05)
table(elite.tmp2$unicox$pvalue < 0.05)
# FALSE  TRUE 
# 641   120
mo.data <- list(omics1 = elite.tmp1$elite.dat,
                omics2 = elite.tmp2$elite.dat)

optk <- getClustNum(data        = mo.data,
                    is.binary   = c(F,F), 
                    try.N.clust = 2:8, # try cluster number from 2 to 8
                    fig.name    = "CLUSTER NUMBER OF TCGA-PRAD")

#the imputed optimal cluster number is 2 arbitrarily

moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("iClusterBayes","SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 2)

save(mo.data,moic.res.list, file = "output_final/1/moic.res.list_2.rda")

load("output_final/1/moic.res.list_2.rda")

cmoic <- getConsensusMOIC(moic.res.list = moic.res.list,
                          fig.name      = "CONSENSUS HEATMAP",
                          distance      = "euclidean",
                          linkage       = "average")

indata <- mo.data
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2), # no truncation for mutation
                     centerFlag = c(T,T), # no center for mutation
                     scaleFlag  = c(T,T)) # no scale for mutation


##
feat   <- moic.res.list$iClusterBayes$feat.res
feat1  <- feat[which(feat$dataset == "omics1"),][1:10,"feature"]
feat2  <- feat[which(feat$dataset == "omics2"),][1:10,"feature"]

annRow <- list(feat1, feat2)
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
col.list   <- list(mRNA.col, lncRNA.col)
##Heatmap
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA"),
             is.binary     = c(F,F), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM"),
             clust.res     = cmoic$clust.res, # consensusMOIC results
             clust.dend    = cmoic$clust.dend, # show no dendrogram for samples
             show.rownames = c(F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             # show.row.dend = c(F,F), # show no dendrogram for features
             annRow        = annRow, # no selected features
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")


###KMplot
surv <- compSurv(moic.res         = cmoic,
                 surv.info        = clin,
                 convt.time       = "m", # convert day unit to month
                 surv.median.line = "h", # draw horizontal line at median survival
                 xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                 fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")
                 
###select more signaficant genes
###Univariate Cox----
load("output_final/1/moic.res.list_2.rda")
mrna <- mo.data$omics1
lnc <- mo.data$omics2
identical(colnames(mrna),colnames(lnc))
a <- rbind(mrna,lnc)

load("all_gene1_data/TCGA.Rda")
clin <- clin[,c(9,10)] %>% na.omit()
identical(colnames(a),rownames(clin))
surv <- clin
surv$BCR=as.numeric(surv$BCR)
surv$BCR.time=as.numeric(surv$BCR.time)

surv.expr <- cbind(surv,t(a))

Coxoutput <- NULL 
library(survival)
for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(BCR.time,BCR) ~ surv.expr[,i], data = surv.expr) 
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}

pcutoff <- 0.01
topgene <- Coxoutput[which(Coxoutput$pvalue < pcutoff),] 
A <- topgene$gene
save(Coxoutput,topgene,file = "output_final/1/gene_70.Rda")
save(A,file = "all_gene1_data/101ML_imput/gene.Rda")
```

## 101ML

```R
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(limma)
library(tidyverse)
library(dplyr)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE) 

# 
load("all_gene1_data/101ML_imput/mm_70.Rda")
mm <- lapply(mm,function(x){
  x[,-c(1:2)] <- scale(x[,-c(1:2)])
  return(x)})

result <- data.frame()
#train dateset
est_data <- mm$TCGA
#valid dateset
val_data_list <- mm
pre_var <- colnames(est_data)[-c(1:2)]

est_dd <- est_data[, c('BCR.time', 'BCR', pre_var)]
val_dd_list <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', pre_var)]})
tune.nodesize(Surv(BCR.time,BCR) ~ .,est_dd)

rf_nodesize <- 5
seed <- 123

##################################
#### 1-1.RSF ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(BCR.time,BCR)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS  = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result, cc)


##################################
#### 1-2.RSF + CoxBoost ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(BCR.time,BCR)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[, 'BCR.time'], est_dd2[, 'BCR'], as.matrix(est_dd2[, -c(1, 2)]), 
                            trace=TRUE, start.penalty = 500, parallel = T)

cv.res <- cv.CoxBoost(est_dd2[, 'BCR.time'], est_dd2[, 'BCR'], as.matrix(est_dd2[, -c(1, 2)]), 
                      maxstepno = 500, K = 10, type = "verweij",  penalty = pen$penalty)
fit <- CoxBoost(est_dd2[, 'BCR.time'], est_dd2[, 'BCR'], as.matrix(est_dd2[, -c(1, 2)]), 
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1],  newstatus = x[, 2], type = "lp")))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CoxBoost')
result <- rbind(result, cc)

##################################
#### 1-3.RSF + Enet ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(BCR.time, BCR)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$BCR.time, est_dd2$BCR))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + ', 'Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

##################################
#### 1-4.RSF + GBM ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(BCR.time, BCR)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize,  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(BCR.time, BCR)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)

# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(BCR.time,BCR)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'GBM')
result <- rbind(result, cc)

##################################
#### 1-5.RSF + Lasso ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(BCR.time, BCR)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$BCR.time, est_dd2$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Lasso')
result <- rbind(result, cc)

##################################
#### 1-6.RSF + plsRcox ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(BCR.time, BCR)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$BCR.time, status = est_dd2$BCR), nt = 10, verbose = FALSE)

fit <- plsRcox(est_dd2[, rid], time = est_dd2$BCR.time, event = est_dd2$BCR, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'plsRcox')
result <- rbind(result, cc)

##################################
#### 1-7.RSF + Ridge ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(BCR.time, BCR)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$BCR.time, est_dd2$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10,
                family = "binomial", alpha = 0)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Ridge')
result <- rbind(result, cc)

##################################
#### 1-8.RSF + StpCox ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(BCR.time,BCR)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(BCR.time, BCR)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}


##################################
#### 1-9.RSF + SuperPC ####
##################################

set.seed(seed)
fit <- rfsrc(Surv(BCR.time,BCR)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$BCR.time,
             censoring.status = est_dd2$BCR, 
             featurenames = colnames(est_dd2)[-c(1, 2)])
set.seed(seed)
fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10, 
                     n.components = 3, 
                     min.features = 2, 
                     max.features = nrow(data$x), 
                     compute.fullcv = TRUE, 
                     compute.preval = TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x = t(w[, -c(1, 2)]), y = w$BCR.time, censoring.status=w$BCR, featurenames = colnames(w)[-c(1, 2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[, 1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'SuperPC')
result <- rbind(result, cc)

##################################
#### 1-10.RSF + survival-SVM ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(BCR.time, BCR)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
fit = survivalsvm(Surv(BCR.time, BCR)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'survival-SVM')
result <- rbind(result,cc)

##################################
#### 2.Enet ####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$BCR.time, est_dd$BCR))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = 'link', newx = as.matrix(x[,-c(1,2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

##################################
#### 3.StepCox ####
##################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(BCR.time,BCR)~., est_dd), direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

for (direction in c("both","backward")) {
  fit <- step(coxph(Surv(BCR.time, BCR)~., est_dd), direction = "both")
  rid <- names(coef(fit))
  est_dd2 <- est_data[,c('BCR.time', 'BCR', rid)]
  val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
  set.seed(seed)
  pen <- optimCoxBoostPenalty(est_dd2[, 'BCR.time'], est_dd2[, 'BCR'], as.matrix(est_dd2[, -c(1,2)]),
                              trace=TRUE, start.penalty = 500, parallel = T)
  cv.res <- cv.CoxBoost(est_dd2[, 'BCR.time'], est_dd2[, 'BCR'], as.matrix(est_dd2[, -c(1,2)]),
                        maxstepno = 500, K = 10 , type = "verweij", penalty = pen$penalty)
  fit <- CoxBoost(est_dd2[, 'BCR.time'], est_dd2[, 'BCR'], as.matrix(est_dd2[, -c(1, 2)]),
                  stepno = cv.res$optimal.step, penalty = pen$penalty)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime=x[, 1], newstatus=x[,2], type="lp")))})
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + CoxBoost')
  result <- rbind(result, cc)
  
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$BCR.time, est_dd2$BCR))
  for (alpha in seq(0.1, 0.9, 0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2, family = "cox",alpha = alpha, nfolds = 10)
    rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
    cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')
    result <- rbind(result, cc)
  }
  set.seed(seed)
  fit <- gbm(formula = Surv(BCR.time, BCR)~., data = est_dd2, distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  # find index for number trees with minimum CV error
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(BCR.time, BCR)~., data = est_dd2, distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 8)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
  cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + GBM')
  result <- rbind(result, cc)
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$BCR.time, est_dd2$BCR))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold=10, 
                  family = "cox", alpha = 1)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Lasso')
  result <- rbind(result, cc)
  set.seed(seed)
  cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$BCR.time, status = est_dd2$BCR), nt = 10, verbose = FALSE)
  fit <- plsRcox(est_dd2[, rid], time = est_dd2$BCR.time,
                 event = est_dd2$BCR, nt = as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1,2)])))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + plsRcox')
  result <- rbind(result, cc)
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$BCR.time, est_dd2$BCR))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold = 10, 
                  family = "cox", alpha = 0)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Ridge')
  result <- rbind(result, cc)
  set.seed(seed)
  fit <- rfsrc(Surv(BCR.time,BCR)~., data = est_dd2,
               ntree = 1000, nodesize = rf_nodesize, 
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + RSF')
  result <- rbind(result, cc)
  data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$BCR.time,
               censoring.status = est_dd2$BCR,
               featurenames = colnames(est_dd2)[-c(1,2)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
  cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                       n.fold = 10,
                       n.components = 3,
                       min.features = 2,
                       max.features = nrow(data$x),
                       compute.fullcv = TRUE,
                       compute.preval = TRUE)
  rs <- lapply(val_dd_list2, function(w){
    test <- list(x = t(w[, -c(1,2)]), y = w$BCR.time, censoring.status = w$BCR, featurenames = colnames(w)[-c(1,2)])
    ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[,1:2], RS = rr)
    return(rr2)
  })
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + SuperPC')
  result <- rbind(result, cc)
  fit = survivalsvm(Surv(BCR.time,BCR)~., data = est_dd2, gamma.mu = 1)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + survival-SVM')
  result <- rbind(result, cc)
}

##################################
#### 4-1.CoxBoost ####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]), maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result, cc)

##################################
#### 4-2.CoxBoost + Enet####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)`!=0), "id"]
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$BCR.time, est_dd2$BCR))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost', ' + Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

##################################
#### 4-3.CoxBoost + GBM####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)`!=0), "id"]
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('BCR.time', 'BCR', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(BCR.time,BCR)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(BCR.time,BCR)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'GBM')
result <- rbind(result, cc)

##################################
#### 4-4.CoxBoost + lasso####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)

cv.res <- cv.CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                stepno=cv.res$optimal.step, penalty=pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$BCR.time, est_dd2$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Lasso')
result <- rbind(result, cc)

##################################
#### 4-5.CoxBoost + plsRcox####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('BCR.time', 'BCR', rid)]})
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$BCR.time, status = est_dd2$BCR), nt = 10, verbose = FALSE)

fit <- plsRcox(est_dd2[, rid], time = est_dd2$BCR.time, event = est_dd2$BCR, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type="lp", newdata = x[, -c(1,2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'plsRcox')
result <- rbind(result, cc)

##################################
#### 4-6.CoxBoost + Ridge####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K=10, type="verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('BCR.time', 'BCR', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$BCR.time, est_dd2$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, 
                family = "cox", alpha = 0)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Ridge')
result <- rbind(result, cc)

##################################
#### 4-7.CoxBoost + StepCox####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('BCR.time', 'BCR', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(BCR.time,BCR)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

##################################
#### 4-8.CoxBoost + SuperPC####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[,-c(1,2)]),
                      maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[,-c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('BCR.time', 'BCR', rid)]})
data <- list(x = t(est_dd2[, -c(1,2)]), y = est_dd2$BCR.time, censoring.status = est_dd2$BCR,
             featurenames = colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 2,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval =TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x=t(w[, -c(1,2)]), y = w$BCR.time, censoring.status = w$BCR, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'SuperPC')
result <- rbind(result, cc)

##################################
#### 4-9.CoxBoost + survival-SVM####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[,-c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
fit = survivalsvm(Surv(BCR.time, BCR)~., data = est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'survival-SVM')
result <- rbind(result, cc)

##################################
#### 4-10.CoxBoost + RSF####
##################################

set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)

cv.res <- cv.CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'BCR.time'], est_dd[, 'BCR'], as.matrix(est_dd[, -c(1,2)]),
                stepno=cv.res$optimal.step, penalty=pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
set.seed(seed)
fit <- rfsrc(Surv(BCR.time,BCR)~., data = est_dd2,
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost', ' + RSF')
result <- rbind(result, cc)
##################################
#### 5.plsRcox####
##################################
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd[,pre_var], time = est_dd$BCR.time, status = est_dd$BCR), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd[,pre_var], time = est_dd$BCR.time, event = est_dd$BCR, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = "lp", newdata = x[, -c(1, 2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result, cc)

##################################
#### 6.superpc####
##################################
data <- list(x = t(est_dd[, -c(1,2)]), y = est_dd$BCR.time, censoring.status = est_dd$BCR, featurenames = colnames(est_dd)[-c(1, 2)])
set.seed(seed) 
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 2,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)
rs <- lapply(val_dd_list, function(w){
  test <- list(x = t(w[,-c(1,2)]), y = w$BCR.time, censoring.status = w$BCR, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result, cc)

##################################
#### 7.GBM####
##################################
set.seed(seed)
fit <- gbm(formula = Surv(BCR.time,BCR)~., data = est_dd, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(BCR.time, BCR)~., data = est_dd, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result, cc)

##################################
#### 8.survival-SVM####
##################################
fit = survivalsvm(Surv(BCR.time,BCR)~., data = est_dd, gamma.mu = 1)
rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('survival - SVM')
result <- rbind(result, cc)

##################################
#### 9.Ridge####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$BCR.time, est_dd$BCR))
set.seed(seed)
fit = glmnet(x1, x2, family = "cox", alpha = 0, lambda = NULL)
cvfit = cv.glmnet(x1, x2,
                  nfold = 10, 
                  family = "cox")

rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = cvfit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time,BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Ridge')
result <- rbind(result, cc)

##################################
####10.Lasso####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$BCR.time, est_dd$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = 'cox', alpha = 1)
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso')
result <- rbind(result, cc)

##################################
#### 10.1.Lasso + CoxBoost####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$BCR.time, est_dd$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1)
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[, 'BCR.time'], est_dd2[, 'BCR'], as.matrix(est_dd2[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd2[, 'BCR.time'], est_dd2[, 'BCR'], as.matrix(est_dd2[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd2[, 'BCR.time'], est_dd2[, 'BCR'], as.matrix(est_dd2[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[,-c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CoxBoost')
result <- rbind(result, cc)

##################################
#### 10.2.Lasso + GBM####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$BCR.time, est_dd$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[,c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('BCR.time', 'BCR', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(BCR.time,BCR)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(BCR.time,BCR)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'GBM')
result <- rbind(result, cc)

##################################
#### 10.3.Lasso + plsRcox####
##################################

x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$BCR.time, est_dd$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1)
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[,c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('BCR.time', 'BCR', rid)]})
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$BCR.time, status = est_dd2$BCR), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd2[, rid], time = est_dd2$BCR.time, event = est_dd2$BCR, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[,-c(1,2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'plsRcox')
result <- rbind(result, cc)

##################################
#### 10.4.Lasso + RSF####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$BCR.time, est_dd$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1)
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid<-rid[-1]
est_dd2 <- est_data[,c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
set.seed(seed)
fit <- rfsrc(Surv(BCR.time,BCR)~., data = est_dd2,
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso', ' + RSF')
result <- rbind(result, cc)

##################################
#### 10.5.Lasso + stepcox####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$BCR.time, est_dd$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1)
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[, c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(BCR.time,BCR)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

##################################
#### 10.6.Lasso + superPC####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$BCR.time, est_dd$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1)
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[,c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('BCR.time', 'BCR', rid)]})
data <- list(x = t(est_dd2[,-c(1,2)]), y = est_dd2$BCR.time, censoring.status = est_dd2$BCR,
             featurenames = colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 2,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x = t(w[,-c(1,2)]), y = w$BCR.time, censoring.status = w$BCR, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[, 1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'SuperPC')
result <- rbind(result, cc)

##################################
#### 10.7.Lasso + survival-SVM####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$BCR.time, est_dd$BCR))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1)
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[,c('BCR.time', 'BCR', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('BCR.time', 'BCR', rid)]})
fit = survivalsvm(Surv(BCR.time,BCR)~., data = est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(BCR.time, BCR) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'survival-SVM')
result <- rbind(result, cc)

save(result,file = "output_final/1/ML/result.Rda")
save(result,file = "output_final/1.2/ML/result.Rda")

load("output_final/1/ML/result.Rda")
result2 <- result

dd2 <- pivot_wider(result2, names_from = 'ID', values_from = 'Cindex') %>% as.data.frame()

dd2[,-1] <- apply(dd2[,-1], 2, as.numeric)

dd2$All <- apply(dd2[,2:5], 1, mean)

dd2$GEO <- apply(dd2[,3:7], 1, mean)

head(dd2)


# write.table(dd2,"outputdata/log2/ML/output_C_index.txt", col.names = T, row.names = F, sep = "\t", quote = F)

# dd2=read.table('output_C_index_39.txt', header=T, sep="\t", check.names=F)


library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

SimpleHeatmap <- function(Cindex_mat, avg_Cindex, 
                          CohortCol, barCol,
                          cellwidth = 1, cellheight = 0.5, 
                          cluster_columns, cluster_rows){
  col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                            col = list("Cohort" = CohortCol),
                            show_annotation_name = F)
  
  row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                            gp = gpar(fill = barCol, col = NA),
                                            add_numbers = T, numbers_offset = unit(-10, "mm"),
                                            axis_param = list("labels_rot" = 0),
                                            numbers_gp = gpar(fontsize = 8, col = "white"),
                                            width = unit(3, "cm")),
                         show_annotation_name = F)
  
  Heatmap(as.matrix(Cindex_mat), name = "AUC",
          right_annotation = row_ha, 
          top_annotation = col_ha,
          col = c("#1CB8B2", "#FFFFFF", "#EEB849"), 
          # col = c("#4195C1", "#FFFFFF", "#CB5746"), 
          rect_gp = gpar(col = "black", lwd = 1),
          cluster_columns = cluster_columns, cluster_rows = cluster_rows, 
          show_column_names = FALSE, 
          show_row_names = TRUE,
          row_names_side = "left",
          width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
          height = unit(cellheight * nrow(Cindex_mat), "cm"),
          column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
          column_title = NULL,
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                      x, y, gp = gpar(fontsize = 8))
          }
  )
}



#dd2 <- dd2[order(dd2$GEO, decreasing = T),]

dt <- dd2[, 3:7]
rownames(dt) <- dd2$Model


Cindex_mat=dt
avg_Cindex <- apply(Cindex_mat, 1, mean)          
avg_Cindex <- sort(avg_Cindex, decreasing = T)     
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]     

avg_Cindex <- as.numeric(format(avg_Cindex, digits = 4, nsmall = 4)) 
if(ncol(Cindex_mat) < 3) { 
  CohortCol <- c("#1CB8B2", "#EEB849")
} else { 
  CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") 
}
names(CohortCol) <- colnames(Cindex_mat)

cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat, 
                    avg_Cindex, 
                    CohortCol, "steelblue",
                    cellwidth = cellwidth, cellheight = cellheight,
                    cluster_columns = F, cluster_rows = F) 

pdf(file.path("output_final/1/ML/Cindex_8.pdf"), width = cellwidth * ncol(Cindex_mat) + 8, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm)
dev.off()
```

## Coefficient and gene signature

```R
load("all_gene1_data/101ML_output/fit.Rda")
a <- fit[["Coeffs"]] %>% as.data.frame()
data <- cbind(rownames(a),a$V1) %>% as.data.frame()
colnames(data) <- c("gene","Coeffs")
data$color <- ifelse(1:23 %in% match(data$gene, c("MGP", "LCN2", "C1QC", "SFRP2")), '#fa9193', '#a1cce2')
data$Coeffs <- as.numeric(data$Coeffs)
data <- data[order(data$Coeffs),]
data$gene <- gsub("_","-",data$gene)

ggplot(data,aes(x=Coeffs,
                y=factor(gene, levels=gene)))+
  theme(legend.title=element_blank())+
  geom_point(aes(size=0.8,
                 color=data$color))+
  geom_segment(aes(x=0,xend=Coeffs,
                   y=gene,yend=gene,
                   color=color,
                   size=0.5))+
  geom_vline(xintercept = 0, linetype="dashed", color="gray",size=0.5) +
  theme_bw(base_rect_size = 1.5)+
  labs(x = "Coeffs", y = NULL, title = "Gene Coeffs") +
  theme(axis.text = element_text(size=10),
        legend.position = 'none',
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=14),
        strip.text = element_text(size=12),
        axis.ticks.x = element_blank())  
```

## ROC

```R
#
bioCol=c("#427AB2",
         "#C59D94",
         "#F09148",
         "#FF9896",
         "#DBDB8D",
         "#AFC7E8")

#
library(pROC)
load("all_gene1_data/101ML_output/RS_107.Rda")
rt <- rs$TCGA
roc=roc(rt$BCR, rt$RS)
aucText=c( paste0("TCGA-AUC=",sprintf("%0.3f",auc(roc))) )
plot(roc,col=bioCol[1],xlim(0,1),ylim(0,1))
for(i in 2:6){
  rt <- rs[[i]]
  roc=roc(rt[,2], rt[,3])
  lines(roc, col=bioCol[i])
  aucText=c(aucText, paste0(names(rs)[i],"-AUC=",sprintf("%0.3f",auc(roc))) )
}
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:6])
dev.off()
```

## compare with other clinical index

```
###TCGA----
load("all_gene1_data/101ML_output/RS_107.Rda")
TCGA <- rs$TCGA
colnames(TCGA) <- c("BCR.time","BCR" , "CRIS")
k <- rownames(TCGA)
load("output_final/0/sgsva_cluster_clin.Rda")
colnames(my)
my <- my[k,-c(1,7)]

identical(rownames(my),rownames(TCGA))
my <- cbind(TCGA,my)
library(compareC)

for (i in 4:ncol(my)) {
  my[,i] <- factor(my[,i])
  
}

tt <- my 
dd <- data.frame()
for (i in colnames(my)[c(3:8)]) {
  fit <- summary(coxph(Surv(BCR.time,BCR)~get(i),tt))
  CC <- fit$concordance[1]%>%as.numeric()
  se <- fit$concordance[2]%>%as.numeric()
  p <- compareC(tt$BCR.time,tt$BCR,tt$CRIS,tt[,i])$pval
  dd <- rbind(dd,data.frame(ID=i,C=CC,SE=se,P=p))
}
ll <- ifelse(dd$P<0.0001,'****',
             ifelse(dd$P<0.001,'***',
                    ifelse(dd$P<0.01,'**',
                           ifelse(dd$P<0.05,'*',''))))
dd$ID
dd$ID <- factor(dd$ID,levels = c("CRIS", "age" , "pT", "pN", "cM",        
                                "gleason_score" ))
ggplot(dd,aes(ID,C,fill=ID))+
  geom_bar(stat='identity',position=position_dodge(0.8),width=0.6)+
  geom_errorbar(aes(ymax=C+1.5*SE,ymin=C-1.5*SE),
                width=0.1,position = position_dodge(0.8),size=0.6)+
  theme_bw(base_rect_size = 1.5)+
  ggtitle("TCGA")+
  labs(x = NULL, y = 'C-index (Compared with IRLS)')+
  theme(axis.text = element_text(size=10),
        legend.position = 'none',
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=14),
        strip.text = element_text(size=12),
        axis.ticks.x = element_blank())+
  scale_fill_npg()+
  scale_y_continuous(expand = c(0,0.01),limits = c(0,1))+
  geom_text(aes(y=0.95,label=ll),size=5)
```

## TCGA-cox

```R
load("all_gene1_data/101ML_output/RS_107.Rda")
TCGA_rs <- rs$TCGA
load("all_gene1_data/TCGA.Rda")
colnames(clin)
clin2 <- clin[,c(1,2,3,4,6)]
identical(rownames(clin2),rownames(TCGA_rs))
colnames(TCGA_rs) <- c("BCR.time","BCR","CRIS")
my <- cbind(TCGA_rs,clin2)
colnames(dat)
my$age<- ifelse(my$age<=65,"age<=65","age>65")
my$age <- factor(my$age,levels = c("age<=65","age>65"))
##针对my进行改造
##设置为因子变量进行多因素回归分析
my$cT <- ifelse(as.numeric(substring(my$cT,2,2))<=2,"<=T2",">T2")
my$cT <- factor(my$cT,levels = c("<=T2",">T2"))
my$pT <- ifelse(as.numeric(substring(my$pT,2,2))<=2,"<=T2",">T2")
my$pT <- factor(my$pT,levels = c("<=T2",">T2"))


Coxoutput <- NULL 
library(survival)
for(i in 3:ncol(my)){
  g <- colnames(my)[i]
  cox <- coxph(Surv(BCR.time,BCR) ~ my[,i], data = my) # 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(id = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           HR.95L = as.numeric(coxSummary$conf.int[,3][1]),
                                           HR.95H = as.numeric(coxSummary$conf.int[,4][1]),
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}

rt <- Coxoutput  
gene <- rt$id
gene
gene[2] <- "age(>65 vs <=65)"
gene[3] <- "pT(T3 vs T2)"
gene[4] <- "pN(N1 vs N0)"
gene[5] <- "cT(T3 vs T2)"


gene
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr," (",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

if(T){
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2))
  
  xlim = c(0,2.6)
  par(mar=c(4,2.5,2,0))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=1
  text(0.5,n:1,gene,adj=0,cex=text.cex)
  text(1.65,n:1,pVal,adj=1,cex=text.cex);
  text(1.5+0.2,n+1,'P-value',cex=text.cex,font=2,adj=1)
  text(2.6,n:1,Hazard.ratio,adj=1,cex=text.cex)
  text(2.45,n+1,'HR (95% CI)',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,0,2,1),mgp=c(2,0.5,0))
  xlim = c(0,10)
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=0,code=3,length=0.05,col="black",lwd=2.5)
  abline(v=1,col="gray",lty=2,lwd=1.5)
  boxcolor = '#67AB9F'
  points(as.numeric(hr), n:1, pch = 15, cex=2,col = boxcolor)
  axis(1)
}




##----
x <- summary(coxph(Surv(BCR.time,BCR)~.,my))
y <- data.frame(id=rownames(x$coefficients),
                HR=x$coefficients[,2],
                HR.95L=x$conf.int[,"lower .95"],
                HR.95H=x$conf.int[,'upper .95'],
                pvalue=x$coefficients[,"Pr(>|z|)"])


y[,-1] <- apply(y[,-1],2,as.numeric)

rt <- y
gene <- rt$id
gene
gene[2] <- "age(>65 vs <=65)"
gene[3] <- "pT(T3 vs T2)"
gene[4] <- "pN(N1 vs N0)"
gene[5] <- "cT(T3 vs T2)"
gene
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr," (",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

if(T){
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2))
  
  xlim = c(0,2.6)
  par(mar=c(4,2.5,2,0))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=1
  text(0.5,n:1,gene,adj=0,cex=text.cex)
  text(1.65,n:1,pVal,adj=1,cex=text.cex);
  text(1.5+0.2,n+1,'P-value',cex=text.cex,font=2,adj=1)
  text(2.6,n:1,Hazard.ratio,adj=1,cex=text.cex)
  text(2.45,n+1,'HR (95% CI)',cex=text.cex,font=2,adj=1,)
  
  par(mar=c(4,0,2,1),mgp=c(2,0.5,0))
  xlim = c(0,6)
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=0,code=3,length=0.05,col="black",lwd=2.5)
  abline(v=1,col="gray",lty=2,lwd=1.5)
  boxcolor = '#67AB9F'
  points(as.numeric(hr), n:1, pch = 15, cex=2,col = boxcolor)
  axis(1)
}
```

## heatmap

```R
library(tidyverse)
load("codes/ML/comPrgML_2.0/Results/cinfo1_final_2-12.Rda")
cinfo1 <- split(cinfo1, cinfo1$Cohort)
dat <- do.call(cbind,lapply(cinfo1, "[",  "pvalue"))
rownames(dat) <- cinfo1$TCGA$method
dat <- as.data.frame(dat)
dat <0.05
dat1 <- ifelse(dat<0.05,"Risky",ifelse(dat > 0.05,"P > 0.05","Protective"))

dat2 <-  as.data.frame(dat1)
colors <- c("#ABD3E1","#FAC795")
library(ggwaffle)
library(reshape)
dat2$id <- rownames(dat2)
dat2$id <- factor(dat2$id,levels = dat2$id)
dat2 <- melt(dat2,id="id")
square_size <- 0.1
ggplot(dat2, aes(id, variable, fill = value)) + 
  geom_waffle(size = square_size)+
  scale_fill_manual(name = "Category",
                    #labels = names(sort_table),
                    values = colors)+
  theme(#panel.border = element_rect(fill=NA,size = 2),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_blank(),
    plot.title = element_text(size = rel(1.2)),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.position = "right")
```

## Signature plot

```R
CohortCol <- brewer.pal(n = length(cinfo), name = "Paired") 
names(CohortCol) <- names(cinfo)

plots <- lapply(cinfo, function(plot.data){
  plot.data$method <- 
    factor(plot.data$method,
           levels = plot.data$method[order(plot.data$C, decreasing = F)])
  
  # compares two concordance indices: the statistical test is a two-sided Student t test for dependent samples.
  C.compare <- plot.data$C[plot.data$method == mySIGname]
  se.compare <- plot.data$se[plot.data$method == mySIGname]
  n.compare <- plot.data$n[plot.data$method == mySIGname]
  RS.compare <- plot.data$RS[plot.data$method == mySIGname][[1]]
  r.combined <- unlist(lapply(plot.data$RS, function(x) cor(x, RS.compare)))
  var.combined <- plot.data$se^2 + se.compare^2 - 2*r.combined*plot.data$se*se.compare
  p <- pt(abs((plot.data$C-C.compare))/(sqrt(var.combined)), n.compare - 1, lower.tail = F) * 2
  plot.data$label <- cut(p, breaks = c(0, 0.05, 0.01, 0.001, 0.0001))
  # plot.data$label <- plyr::mapvalues(x = plot.data$label,
  #                                    from = c("(0,0.0001]", "(0.0001,0.001]", "(0.001,0.01]", "(0.01,0.05]"), 
  #                                    to = c("****", "***", "**", "*"))
  
  return(ggplot(plot.data, aes(x = method, y = C, fill = Cohort)) +
           geom_errorbar(aes(ymin = C - 1.96 * se, ymax = C + 1.96 * se), width = .1) +
           geom_point(color = CohortCol[unique(plot.data$Cohort)], size = 2.5) +
           # geom_text(aes(x = method, y = max(plot.data$C + 1.96 * plot.data$se - 0.05), label = label)) +
           geom_hline(yintercept = 0.6, linetype = "dashed") +
           ggtitle(label = unique(plot.data$Cohort)) +
           coord_flip() + 
           theme_classic() +
           theme(panel.border = element_rect(fill = NA, size = 1),
                 axis.title = element_blank(),
                 legend.position = "none"))
})

plot_grid(plotlist = plots, nrow = 1)
ggsave(file.path(fig.path, "comparison_final_2.pdf"), width = 20, height = 15)
```

