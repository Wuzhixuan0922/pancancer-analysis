rm(list = ls())
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(clusterProfiler)
library(pheatmap)
library(tidyverse)
library(survminer)
library(survival)
library(limma)
library(ggpubr)
library(stringr)
library(RColorBrewer)
library(export)
library(survivalROC)
library(pheatmap)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(ReactomePA)
library(igraph)
library(ggraph)
library(forestplot)
library(ggradar)
library(fmsb)
library(circlize)
library(ggsci)
library(parallel)
library(maftools)
library(patchwork)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = "en") 
gene="XXX"  #需要替换



VVV<-c( 'ACC', 'BLCA' ,'BRCA' ,'CESC', 'CHOL','COAD','READ','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG',
        'LIHC','LUAD', 'LUSC','OV','PAAD','PCPG', 'PRAD','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS' )

names<-read.table(paste0(gene,".txt"),header = F)
names<-names$V1

# aaa<-dplyr::filter(survland.coxplot1,sum>=4)
# names<-aaa$name
# names<-as.character(names)

heatmapdata<-read.csv(paste0(VVV[1],"差异分析结果.csv"),header = T,row.names = 1)
names1<-intersect(names,rownames(heatmapdata))

length(names1)==length(names)

for (i in 1:31) {
  print(i)
  heatmapdata<-read.csv(paste0(VVV[i],"差异分析结果.csv"),header = T,row.names = 1)
  # rownames(heatmapdata)[which(rownames(heatmapdata)=="C9orf41")]<-"CARNMT1"
  # rownames(heatmapdata)[which(rownames(heatmapdata)=="LINC00116")]<-"MTLN"
  # rownames(heatmapdata)[which(rownames(heatmapdata)=="TMEM173")]<-"STING1"
  # rownames(heatmapdata)[which(rownames(heatmapdata)=="H2AFX")]<-"H2AX"
  # rownames(heatmapdata)[which(rownames(heatmapdata)=="ERO1L")]<-"ERO1A"
  # rownames(heatmapdata)[which(rownames(heatmapdata)=="IARS")]<-"IARS1"
  heatmapdata<-heatmapdata[names,]
  heatmapdata$ID<-rownames(heatmapdata)
  aaaaa<-dplyr::select(heatmapdata,ID,logFC)
  bbbbb<-dplyr::select(heatmapdata,ID,P.Value)
  colnames(aaaaa)[2]<-VVV[i]
  colnames(bbbbb)[2]<-VVV[i]
  lofFCdata<-inner_join(aaaaa,lofFCdata,by="ID")
  pvaluedata<-inner_join(bbbbb,pvaluedata,by="ID")
}
rownames(lofFCdata)<-lofFCdata$ID
lofFCdata<-lofFCdata[,-1]
rownames(pvaluedata)<-pvaluedata$ID
pvaluedata<-pvaluedata[,-1]

#为确保热图基因的顺序和分组文件的顺序是一致的
#对logfc和p.val两个文件排序
logfcOrdered <- lofFCdata[names,]
p.valOrdered <- pvaluedata[names,]

## 对logfc进行分类
logfcCat <- apply(logfcOrdered, 2, function(x){
  cut(x, breaks = c(-Inf, -2, -1, 0,1, 2, Inf),
      labels = c("< -2", "-2 - -1", "-1 - 0", "0 - 1", "1 - 2", "> 2"))
})
rownames(logfcCat) <- rownames(logfcOrdered)

## 确保两个数据集的列名和行名顺序一致
p.valOrdered<-p.valOrdered[,VVV]
p.valOrdered<-p.valOrdered[names,]
logfcCat<-logfcCat[,VVV]
logfcCat<-logfcCat[names,]

logfcCat[p.valOrdered >= 0.05] <- "P >= 0.05"


unique(matrix(logfcCat, ncol  = 1))


col_cat <- c("> 2" = "#A80C3A", "1 - 2" = "#ED5E57", "0 - 1" = "#fab2a4","-1 - 0" = "#DDD3D2", 
             "-2 - -1" = "#6B9AB7", "< -2" = "#2F5B89", "P >= 0.05" = "white")
cell_fun <- function(logfc, dataP, logfcCutoff = 1, PCutoff = 0.05, 
                     darkcol = "black", lightcol = "white", digit = 2, fontsize  = 6){
  function(j, i, x, y, width, height, fill){
    if(abs(logfc[i,j]) > logfcCutoff & dataP[i,j] < PCutoff){
      grid.text(round(logfc, digit)[i, j], x, y, 
                gp = gpar(fontsize = fontsize, col  = lightcol))
    }else{
      grid.text(round(logfc, digit)[i, j], x, y, 
                gp = gpar(fontsize = fontsize, col  = darkcol))
    }
  }
}

logfcCat1<-lofFCdata
logfcCat1$sum<-rowSums(logfcCat1)
logfcCat1<-arrange(logfcCat1,desc(sum))

logfcCat<-logfcCat[rownames(logfcCat1),]
logfcOrdered<-logfcOrdered[rownames(logfcCat1),]
p.valOrdered<-p.valOrdered[rownames(logfcCat1),]
logfcOrdered<-logfcOrdered[,VVV]
p.valOrdered<-p.valOrdered[,VVV]

pdf("rankHeatmap.pdf", width = 10, height = 4)
Heatmap(matrix = logfcCat, 
        name = "logFC", #主要图例的标题
        rect_gp = gpar(col = "grey", lwd = 1), #不画边框，或者用col = "grey"画灰色边框
        col = col_cat, #热图颜色
        row_names_side = "left", 
        cell_fun = cell_fun(logfcOrdered, p.valOrdered), 
        row_names_gp = gpar(fontsize = 8), #基因名字号
        column_names_gp = gpar(fontsize = 8), #肿瘤类型字号
        column_names_rot = 90)

dev.off()


geneSets <- getGmt(paste0(gene,".gmt"))
names<-read.table(paste0(gene,".txt"),header = F)
names<-names$V1

load("pancancer_drawdata.RDATA")
AAA<-c("ID",colnames(drawdata)[8:21097])

# load("pancancer_drawdata.RDATA")
# AAA<-c("ID",colnames(drawdata)[8:19616])

drawdata<-drawdata[,AAA]

names[which(names %in% colnames(drawdata)==FALSE)]

rownames(drawdata)<-drawdata$ID
drawdata<-drawdata[,-1]
drawdata<-as.data.frame(t(drawdata))
drawdata<-as.matrix(drawdata)
ssgsea_results <- GSVA::gsva(drawdata,geneSets,method="ssgsea",verbose=T,ssgsea.norm=TRUE)
res<-as.data.frame(t(ssgsea_results))
res1<-res

write.csv(res1,paste0(gene,"_ssgsea.csv"))
save(res1,file =paste0(gene,"_ssgsea.RDATA") )

#load(paste0(gene,"_ssgsea.RDATA"))
#res=res1
res$ID<-row.names(res)
row.names(res)<-NULL
colnames(res)[1]<-gene
res<-dplyr::select(res,"ID",gene)


load("pancancer_drawdata.Rdata")

drawdata<-inner_join(res,drawdata,by='ID')
drawdata<-dplyr::select(drawdata,ID,Cancer,gender,Stage,status,time,Type,gene,everything())

aaa<-dplyr::select(drawdata,ID,gene)

load('Stagedataa.RDATA')
load('OSdata.RDATA')
load('immunedata.RDATA')
load('MSIDATA.RDATA')
load("TMEscores.Rdata")
load("estimateScoresdata.RDATA")


Stagedataa<-inner_join(Stagedataa,aaa, by = "ID")
colnames(aaa)[1]<-"ID.1"
exprSet_pair<-inner_join(exprSet_pair,aaa, by = "ID.1")
colnames(aaa)[1]<-"ID"
OSdata<-drawdata

immunedata<-inner_join(immunedata,aaa, by = "ID")
MSIDATA<-inner_join(MSIDATA,aaa, by = "ID")
TMEscores<-inner_join(TMEscores,aaa, by = "ID")
estimateScoresdata<-inner_join(estimateScoresdata,aaa, by = "ID")


load('DFI_DATA.RDATA')
load('DSS_DATA.RDATA')
load('PFI_DATA.RDATA')
load('TMBDATA.RDATA')
load("ImmuneCelldata.RDATA")
aaa$ID<-str_replace_all(aaa$ID,"-",".")
DFI_DATA<-inner_join(DFI_DATA,aaa, by = "ID")
DSS_DATA<-inner_join(DSS_DATA,aaa, by = "ID")
PFI_DATA<-inner_join(PFI_DATA,aaa, by = "ID")
TMBDATA<-inner_join(TMBDATA,aaa, by = "ID")
ImmuneCelldata<-inner_join(ImmuneCelldata,aaa, by = "ID")


drawdata111<-drawdata
drawdata111$ID<-str_replace_all(drawdata111$ID,  '-' ,'.'  )
GSVAResult<-read.csv("hallmark_ssgsea_all.csv",header = T,sep = ",")
rownames(GSVAResult)<-GSVAResult$ID
GSVAResult<-GSVAResult[,-1]
GSVAResult<-as.data.frame(t(GSVAResult))
Index<-c( 'ACC', 'BLCA' ,'BRCA' ,'CESC', 'CHOL','COAD','READ','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG',
          'LIHC','LUAD', 'LUSC','MESO','OV','PAAD','PCPG', 'PRAD','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM'  )



index<-vector()
name<-vector()
for (i in 1:33) {
  aaa<-filter(drawdata, Cancer == unique(drawdata$Cancer)[i])
  index<-c(mean(aaa[,gene]),index)
  name<-c(unique(drawdata$Cancer)[i],name)
}
dat<-data.frame(name,index)
value<-as.character(arrange(dat,index)$name)
aaaa<-drawdata
aaaa$Cancer<-factor(aaaa$Cancer,levels =value )
p1<-ggplot(data=aaaa,aes(x=Cancer,y=get(gene),fill=Cancer))+
  geom_boxplot()+
  ylab(label =paste0(gene," score") )+
  theme(legend.position = 'none' ,axis.text.x =element_text(angle=90,size=8,hjust =1,vjust =0.5,colour = "black" ),
        axis.text.y =element_text(colour = "black" ))

ggsave(p1,filename = paste0(gene,"_泛癌表达.pdf"),height = 4,width = 8)
ggsave(p1,filename = paste0(gene,"_泛癌表达.tiff"),height = 4,width = 8)
  

Stagedataa11<-Stagedataa
Stagedataa11$Stage[Stagedataa11$Stage=="Stage I"]<-"Stage I+II"
Stagedataa11$Stage[Stagedataa11$Stage=="Stage II"]<-"Stage I+II"
Stagedataa11$Stage[Stagedataa11$Stage=="Stage III"]<-"Stage III+IV"
Stagedataa11$Stage[Stagedataa11$Stage=="Stage IV"]<-"Stage III+IV"

drawdata$Type<-factor(drawdata$Type,levels = c('Tumor','Normal'))
my_comparisons <- list( c("Tumor", "Normal") )

Stagedataa$Stage<-factor(Stagedataa$Stage,levels=c('Stage I','Stage II',"Stage III","Stage IV"))
my_comparisons1 <- list( c("Stage I", "Stage II"),c("Stage I", "Stage III"),c("Stage I", "Stage IV"),c("Stage II", "Stage III"),c("Stage II", "Stage IV"),c("Stage III", "Stage IV") )

exprSet_pair$Type<-factor(exprSet_pair$Type,levels = c('Tumor','Normal'))
my_comparisons2 <- list( c("Tumor", "Normal") )

Stagedataa11$Stage<-factor(Stagedataa11$Stage,levels=c('Stage I+II','Stage III+IV'))
my_comparisons3 <- list( c("Stage I+II", "Stage III+IV"))

for (a in 1:33) {
  print(a)
  p <- ggboxplot(filter(drawdata,Cancer==Index[a]), x = "Type", y = gene,
                 fill = "Type",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
    theme(legend.position='none')+
    ylab(label = paste0(gene," score") )+
    xlab(label = paste(Index[a]))+
    stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
  ggsave(p,filename = paste(gene,Index[a],'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)
  p <- ggboxplot(filter(Stagedataa, Cancer==Index[a],Type == "Tumor"), x = "Stage", y = gene,
                 fill = "Stage",legend=F,bxp.errorbar=T)+
    theme(legend.position='none')+
    ylab(label = paste0(gene," score") )+
    xlab(label = paste(Index[a]))+
    stat_compare_means(comparisons = my_comparisons1,method = "t.test",aes(label = ..p.signif..))
  ggsave(p,filename = paste(gene,Index[a],'分期差异表达.pdf',sep=' '),width = 5.6,height = 4.22)
  p <- ggboxplot(filter(Stagedataa, Cancer==Index[a],Type == "Tumor"), x = "Stage", y = gene,
                 fill = "Stage",legend=F,bxp.errorbar=T)+
    theme(legend.position='none')+
    ylab(label = paste0(gene," score") )+
    xlab(label = paste(Index[a]))+
    stat_compare_means(comparisons = my_comparisons3,method = "t.test",aes(label = ..p.signif..))
  ggsave(p,filename = paste(gene,Index[a],'分期差异表达.pdf',sep=' '),width = 5.6,height = 4.22)
  p<-ggpaired(filter(exprSet_pair,Cancer==Index[a]), x = 'Type', y = gene,id='ID',
              color = 'black',fill = "Type",  palette =c("#E7B800", "#00AFBB"), line.color = "gray", line.size = 0.4,short.panel.labs = FALSE)+
    stat_compare_means(aes(label = ..p.signif..),method = "t.test",paired = TRUE, comparisons = my_comparisons2)+
    ylab(label = paste0(gene," score") )+
    xlab(label = paste(Index[a]))

  indexx<-filter(OSdata,Cancer==Index[a])
  med.exp<-median(indexx[,gene])
  more.med.exp.index<-which(indexx[,gene]>=med.exp)
  less.med.exp.index<-which(indexx[,gene]<med.exp)
  indexx[,gene][more.med.exp.index]<-paste0('High expression (',
                                            length(more.med.exp.index),')')
  indexx[,gene][less.med.exp.index]<-paste0('Low expression (',
                                            length(less.med.exp.index),')')
  Gene<-indexx[,gene]
  indexx.fit<-survfit(Surv(time, status) ~ Gene, data = indexx)
  p<-ggsurvplot(indexx.fit,
                data=indexx,
                risk.table = T,
                risk.table.col = 'strata',
                risk.table.y.text = F,
                risk.table.title = 'Number at risk',
                pval = TRUE,
                #conf.int = TRUE,
                xlab = 'Time (days)',
                ggtheme = theme_light(),
                surv.median.line = 'hv',
                title=paste(Index[a], gene, 'Survival',sep = ' ' ))
  pdf(file = paste( Index[[a]],gene,'生存分析.pdf',sep=' '),width = 5.6,height = 5.34,onefile=F )
  print(p)
  dev.off()
}


for (i in 1:33) {
  print(i)
  test<-dplyr::filter(OSdata,Cancer == Index[[i]])
  res.cut <- surv_cutpoint(test, time = "time", event = "status",
                           variables = gene)
  summary(res.cut)
  #Categorize variables
  res.cat <- surv_categorize(res.cut)
  head(res.cat)
  #Fit survival curves and visualize
  fit <- survfit(Surv(time, status) ~get(gene), data = res.cat)
  p<-ggsurvplot(fit,
                data=res.cat,
                risk.table = T,
                risk.table.col = 'strata',
                risk.table.y.text = F,
                risk.table.title = 'Number at risk',
                pval = TRUE,
                #conf.int = TRUE,
                xlab = 'Time (days)',
                ggtheme = theme_light(),
                surv.median.line = 'hv',
                title=paste( Index[[i]], gene, 'Survival',sep = ' ' ))
  pdf(file = paste( Index[[i]],gene,'最佳cutoff生存分析.pdf',sep=' '),width = 5.6,height = 5.34,onefile=F )
  print(p)
  dev.off()
}


for (i in 1:33) {
  if(Index[[i]] %in%  unique(DFI_DATA$Cancer)  &  Index[[i]] != 'GBM') {#GBM只有1个DFI病例，做生存分析会出错
    indexx<-filter(DFI_DATA,Cancer==Index[[i]] )
    med.exp<-median(indexx[,gene])
    more.med.exp.index<-which(indexx[,gene]>=med.exp)
    less.med.exp.index<-which(indexx[,gene]<med.exp)
    indexx[,gene][more.med.exp.index]<-paste0('High expression (',
                                              length(more.med.exp.index),')')
    indexx[,gene][less.med.exp.index]<-paste0('Low expression (',
                                              length(less.med.exp.index),')')
    Gene<-indexx[,gene]
    indexx.fit<-survfit(Surv(DFI.time, DFI) ~ Gene, data = indexx)
    p<-ggsurvplot(indexx.fit,
                  data=indexx,
                  risk.table = T,
                  risk.table.col = 'strata',
                  risk.table.y.text = F,
                  risk.table.title = 'Number at risk',
                  pval = TRUE,
                  #conf.int = TRUE,
                  xlab = 'Time (days)',
                  ggtheme = theme_light(),
                  surv.median.line = 'hv',
                  title=paste(as.character(Index[[i]]) , gene, 'DFI Survival',sep = ' ' ))
    pdf(file = paste( Index[[i]],gene,'DFI生存分析.pdf',sep=' '),width = 5.6,height = 5.34,onefile=F )
    print(p)
    dev.off()
  }
  print(i)
}

for (i in 1:33) {
  if(Index[[i]] %in%  unique(DFI_DATA$Cancer)  &  Index[[i]] != 'GBM') {#GBM只有1个DFI病例，做生存分析会出错
    test<-dplyr::filter(DFI_DATA,Cancer == Index[[i]])
    res.cut <- surv_cutpoint(test, time = "DFI.time", event = "DFI",
                             variables = gene)
    summary(res.cut)
    #Categorize variables
    res.cat <- surv_categorize(res.cut)
    head(res.cat)
    #Fit survival curves and visualize
    fit <- survfit(Surv(DFI.time, DFI) ~get(gene), data = res.cat)
    p<-ggsurvplot(fit,
                  data=res.cat,
                  risk.table = T,
                  risk.table.col = 'strata',
                  risk.table.y.text = F,
                  risk.table.title = 'Number at risk',
                  pval = TRUE,
                  #conf.int = TRUE,
                  xlab = 'Time (days)',
                  ggtheme = theme_light(),
                  surv.median.line = 'hv',
                  title=paste( Index[[i]], gene, 'Survival',sep = ' ' ))
    pdf(file = paste( Index[[i]],gene,'最佳cutoff生存分析.pdf',sep=' '),width = 5.6,height = 5.34,onefile=F )
    print(p)
    dev.off()
  }
  print(i)
}

#------------------------------------------------------------------
for (i in 1:33) {
  if(Index[[i]] %in%  unique(DSS_DATA$Cancer) ) {
    indexx<-filter(DSS_DATA,Cancer==Index[[i]] )
    med.exp<-median(indexx[,gene])
    more.med.exp.index<-which(indexx[,gene]>=med.exp)
    less.med.exp.index<-which(indexx[,gene]<med.exp)
    indexx[,gene][more.med.exp.index]<-paste0('High expression (',
                                              length(more.med.exp.index),')')
    indexx[,gene][less.med.exp.index]<-paste0('Low expression (',
                                              length(less.med.exp.index),')')
    Gene<-indexx[,gene]
    indexx.fit<-survfit(Surv(DSS.time, DSS) ~ Gene, data = indexx)
    p<-ggsurvplot(indexx.fit,
                  data=indexx,
                  risk.table = T,
                  risk.table.col = 'strata',
                  risk.table.y.text = F,
                  risk.table.title = 'Number at risk',
                  pval = TRUE,
                  #conf.int = TRUE,
                  xlab = 'Time (days)',
                  ggtheme = theme_light(),
                  surv.median.line = 'hv',
                  title=paste(as.character(Index[[i]]) , gene, 'DSS Survival',sep = ' ' ))
    pdf(file = paste( Index[[i]],gene,'DSS生存分析.pdf',sep=' '),width = 5.6,height = 5.34,onefile=F )
    print(p)
    dev.off()
  }
  print(i)
}

for (i in 1:33) {
  if(Index[[i]] %in%  unique(DSS_DATA$Cancer)) {
    test<-dplyr::filter(DSS_DATA,Cancer == Index[[i]])
    res.cut <- surv_cutpoint(test, time = "DSS.time", event = "DSS",
                             variables = gene)
    summary(res.cut)
    #Categorize variables
    res.cat <- surv_categorize(res.cut)
    head(res.cat)
    #Fit survival curves and visualize
    fit <- survfit(Surv(DSS.time, DSS) ~get(gene), data = res.cat)
    p<-ggsurvplot(fit,
                  data=res.cat,
                  risk.table = T,
                  risk.table.col = 'strata',
                  risk.table.y.text = F,
                  risk.table.title = 'Number at risk',
                  pval = TRUE,
                  #conf.int = TRUE,
                  xlab = 'Time (days)',
                  ggtheme = theme_light(),
                  surv.median.line = 'hv',
                  title=paste( Index[[i]], gene, 'Survival',sep = ' ' ))
    pdf(file = paste( Index[[i]],gene,'最佳cutoff生存分析.pdf',sep=' '),width = 5.6,height = 5.34,onefile=F )
    print(p)
    dev.off()
  }
  print(i)
}

#------------------------------------------------------------------
for (i in 1:33) {
  if(Index[[i]] %in%  unique(PFI_DATA$Cancer) ) {
    indexx<-filter(PFI_DATA,Cancer==Index[[i]] )
    med.exp<-median(indexx[,gene])
    more.med.exp.index<-which(indexx[,gene]>=med.exp)
    less.med.exp.index<-which(indexx[,gene]<med.exp)
    indexx[,gene][more.med.exp.index]<-paste0('High expression (',
                                              length(more.med.exp.index),')')
    indexx[,gene][less.med.exp.index]<-paste0('Low expression (',
                                              length(less.med.exp.index),')')
    Gene<-indexx[,gene]
    indexx.fit<-survfit(Surv(PFI.time, PFI) ~ Gene, data = indexx)
    p<-ggsurvplot(indexx.fit,
                  data=indexx,
                  risk.table = T,
                  risk.table.col = 'strata',
                  risk.table.y.text = F,
                  risk.table.title = 'Number at risk',
                  pval = TRUE,
                  #conf.int = TRUE,
                  xlab = 'Time (days)',
                  ggtheme = theme_light(),
                  surv.median.line = 'hv',
                  title=paste(as.character(Index[[i]]) , gene, 'PFI Survival',sep = ' ' ))
    pdf(file = paste( Index[[i]],gene,'PFI生存分析.pdf',sep=' '),width = 5.6,height = 5.34,onefile=F )
    print(p)
    dev.off()
  }
  print(i)
}

for (i in 1:33) {
  if(Index[[i]] %in%  unique(PFI_DATA$Cancer)) {
    test<-dplyr::filter(PFI_DATA,Cancer == Index[[i]])
    res.cut <- surv_cutpoint(test, time = "PFI.time", event = "PFI",
                             variables = gene)
    summary(res.cut)
    #Categorize variables
    res.cat <- surv_categorize(res.cut)
    head(res.cat)
    #Fit survival curves and visualize
    fit <- survfit(Surv(PFI.time, PFI) ~get(gene), data = res.cat)
    p<-ggsurvplot(fit,
                  data=res.cat,
                  risk.table = T,
                  risk.table.col = 'strata',
                  risk.table.y.text = F,
                  risk.table.title = 'Number at risk',
                  pval = TRUE,
                  #conf.int = TRUE,
                  xlab = 'Time (days)',
                  ggtheme = theme_light(),
                  surv.median.line = 'hv',
                  title=paste( Index[[i]], gene, 'Survival',sep = ' ' ))
    pdf(file = paste( Index[[i]],gene,'最佳cutoff生存分析.pdf',sep=' '),width = 5.6,height = 5.34,onefile=F )
    print(p)
    dev.off()
  }
  print(i)
}



#--------------------------------------------------------------------------------------------------------------------------------

value<-c(colnames(immunedata)[19619:19644])
genecor.parallel <- function(data,gene,cl){
  cl <- makeCluster(cl)
  y <- as.numeric(data[gene,])
  rownames <- rownames(data)
  dataframes <- do.call(rbind, parLapply(cl=cl,rownames, function(x){
    dd  <- cor.test(as.numeric(data[x,]), y, type="spearman")
    data.frame(Gene_1=gene, Gene_2=x, cor=dd$estimate, p.value=dd$p.value)
  }))
  stopCluster(cl)
  return(dataframes)
}
genecor_circleplot <- function(x){
  Corr <- data.frame(rbind(data.frame(Gene=x[,1], Correlation=x[,3]), 
                           data.frame(Gene=x[,2], Correlation=x[,3])), stringsAsFactors = F)      
  Corr$Index <- seq(1,nrow(Corr),1) #记录基因的原始排序，记录到Index列
  Corr <- Corr[order(Corr[,1]),] #按照基因名排序
  corrsp <- split(Corr,Corr$Gene)
  corrspe <- lapply(corrsp, function(x){x$Gene_Start<-0
  
  if (nrow(x)==1){x$Gene_End<-1}else{
    x$Gene_End<-sum(abs(x$Correlation))} 
  x})
  GeneID <- do.call(rbind,corrspe)
  GeneID <- GeneID[!duplicated(GeneID$Gene),]
 
  mycol <- c(pal_d3("category20c")(20),pal_d3("category20")(20))
  n <- nrow(GeneID)
  GeneID$Color <- mycol[1:n]
  
  Corr[,2] <- abs(Corr[,2]) 
  corrsl <- split(Corr,Corr$Gene)
  aaaaa <- c()
  corrspl <- lapply(corrsl,function(x){nn<-nrow(x)
  for (i in 1:nn){
    aaaaa[1] <- 0
    aaaaa[i+1] <- x$Correlation[i]+aaaaa[i]}
  bbbbb <- data.frame(V4=aaaaa[1:nn],V5=aaaaa[2:(nn+1)])
  bbbbbb <- cbind(x,bbbbb)
  bbbbbb
  })
  Corr <- do.call(rbind,corrspl)
  
  Corr <- Corr[order(Corr$Index),]

  #把它写入Links里，start_1和end_1对应Gene_1，start_2和end_2对应Gene_2
  x$start_1 <- Corr$V4[1:(nrow(Corr)/2)]
  x$end_1 <- Corr$V5[1:(nrow(Corr)/2)]
  x$start_2 <- Corr$V4[(nrow(Corr)/2 + 1):nrow(Corr)]
  x$end_2 <- Corr$V5[(nrow(Corr)/2 + 1):nrow(Corr)]
  
  color <- data.frame(colorRampPalette(c("#67BE54", "#FFFFFF", "#F82C2B"))(201))
  for (i in 1:nrow(x)){
    x[i,8] <- substring(color[x[i,3] * 100 + 101, 1], 1, 7)
  }
  names(x)[8] <- "color"

  #par(mar=rep(0,4))
  circos.clear()
  circos.par(start.degree = 90, #从哪里开始画，沿着逆时针顺序
             gap.degree = 5, #基因bar之间的间隔大小
             track.margin = c(0,0.23), #值越大，基因跟连线的间隔越小
             cell.padding = c(0,0,0,0)
  )
  circos.initialize(factors = GeneID$Gene,
                    xlim = cbind(GeneID$Gene_Start, GeneID$Gene_End))
  
  #先画基因
  circos.trackPlotRegion(ylim = c(0, 1), factors = GeneID$Gene, 
                         track.height = 0.05, #基因线条的胖瘦
                         panel.fun = function(x, y) {
                           name = get.cell.meta.data("sector.index") 
                           i = get.cell.meta.data("sector.numeric.index") 
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           circos.text(x = mean(xlim), y = 1,
                                       labels = name,
                                       cex = 0.6, #基因ID文字大小
                                       niceFacing = TRUE, #保持基因名的头朝上
                                       facing = "bending", #基因名沿着圆弧方向，还可以是reverse.clockwise
                                       adj = c(0.5, -2.8), #基因名所在位置，分别控制左右和上下
                                       font = 2 #加粗
                           )
                           circos.rect(xleft = xlim[1], 
                                       ybottom = ylim[1],
                                       xright = xlim[2], 
                                       ytop = ylim[2],
                                       col = GeneID$Color[i],
                                       border = GeneID$Color[i])
                           
                           circos.axis(labels.cex = 0.7, 
                                       direction = "outside"
                           )})
  
  #画连线
  for(i in 1:nrow(x)){
    circos.link(sector.index1 = x$Gene_1[i], 
                point1 = c(x[i, 4], x[i, 5]),
                sector.index2 = x$Gene_2[i], 
                point2 = c(x[i, 6], x[i, 7]),
                col = paste(x$color[i], "C9", sep = ""), 
                border = FALSE, 
                rou = 0.7
    )}
  
  #画图例
  i <- seq(0,0.995,0.005)
  rect(-1+i/2, #xleft
       -1, #ybottom
       -0.9975+i/2, #xright
       -0.96, #ytop
       col = paste(as.character(color[,1]), "FF", sep = ""),
       border = paste(as.character(color[,1]), "FF", sep = ""))
  text(-0.97, -1.03, "-1")
  text(-0.51, -1.03, "1")
  
}
# 画图的函数
genecor_circleplot1 <- function(x){
  Corr <- data.frame(rbind(data.frame(Gene=x[,1], Correlation=x[,3]), 
                           data.frame(Gene=x[,2], Correlation=x[,3])), stringsAsFactors = F)      
  Corr$Index <- seq(1,nrow(Corr),1) #记录基因的原始排序，记录到Index列
  Corr <- Corr[order(Corr[,1]),] #按照基因名排序
  corrsp <- split(Corr,Corr$Gene)
  corrspe <- lapply(corrsp, function(x){x$Gene_Start<-0
  
  #依次计算每个基因的相关系数总和，作为基因终止位点
  if (nrow(x)==1){x$Gene_End<-1}else{
    x$Gene_End<-sum(abs(x$Correlation))} 
  x})
  GeneID <- do.call(rbind,corrspe)
  GeneID <- GeneID[!duplicated(GeneID$Gene),]
  
  #基因配色
  mycol <- c(pal_d3("category20c")(20),pal_d3("category20")(20))
  n <- nrow(GeneID)
  GeneID$Color <- mycol[1:n]
  
  #连线的宽度是相关系数的绝对值
  Corr[,2] <- abs(Corr[,2]) 
  corrsl <- split(Corr,Corr$Gene)
  aaaaa <- c()
  corrspl <- lapply(corrsl,function(x){nn<-nrow(x)
  for (i in 1:nn){
    aaaaa[1] <- 0
    aaaaa[i+1] <- x$Correlation[i]+aaaaa[i]}
  bbbbb <- data.frame(V4=aaaaa[1:nn],V5=aaaaa[2:(nn+1)])
  bbbbbb <- cbind(x,bbbbb)
  bbbbbb
  })
  Corr <- do.call(rbind,corrspl)
  
  #根据Index列，把基因恢复到原始排序
  Corr <- Corr[order(Corr$Index),]
  
  #V4是起始位置，V5是终止位置
  #把它写入Links里，start_1和end_1对应Gene_1，start_2和end_2对应Gene_2
  x$start_1 <- Corr$V4[1:(nrow(Corr)/2)]
  x$end_1 <- Corr$V5[1:(nrow(Corr)/2)]
  x$start_2 <- Corr$V4[(nrow(Corr)/2 + 1):nrow(Corr)]
  x$end_2 <- Corr$V5[(nrow(Corr)/2 + 1):nrow(Corr)]
  
  #连线（相关系数）的配色
  #相关系数最大为1，最小-1，此处设置201个颜色
  #-1到0就是前100，0到1就是后100
  color <- data.frame(colorRampPalette(c("#67BE54", "#FFFFFF", "#F82C2B"))(201))
  #根据相关系数的数值，给出相应的颜色
  for (i in 1:nrow(x)){
    x[i,8] <- substring(color[x[i,3] * 100 + 101, 1], 1, 7)
  }
  names(x)[8] <- "color"
  
  #绘图区设置
  #par(mar=rep(0,4))
  circos.clear()
  circos.par(start.degree = 90, #从哪里开始画，沿着逆时针顺序
             gap.degree = 5, #基因bar之间的间隔大小
             track.margin = c(0,0.23), #值越大，基因跟连线的间隔越小
             cell.padding = c(0,0,0,0)
  )
  circos.initialize(factors = GeneID$Gene,
                    xlim = cbind(GeneID$Gene_Start, GeneID$Gene_End))
  
  #先画基因
  circos.trackPlotRegion(ylim = c(0, 1), factors = GeneID$Gene, 
                         track.height = 0.05, #基因线条的胖瘦
                         panel.fun = function(x, y) {
                           name = get.cell.meta.data("sector.index") 
                           i = get.cell.meta.data("sector.numeric.index") 
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           circos.text(x = mean(xlim), y = 1,
                                       labels = name,
                                       cex = 0.4, #基因ID文字大小
                                       niceFacing = TRUE, #保持基因名的头朝上
                                       facing = "reverse.clockwise", #基因名沿着圆弧方向，还可以是reverse.clockwise
                                       adj = c(1.2,0), #基因名所在位置，分别控制左右和上下
                                       font = 2 #加粗
                           )
                           circos.rect(xleft = xlim[1], 
                                       ybottom = ylim[1],
                                       xright = xlim[2], 
                                       ytop = ylim[2],
                                       col = GeneID$Color[i],
                                       border = GeneID$Color[i])
                           
                           circos.axis(labels.cex = 0.7, 
                                       direction = "outside"
                           )})
  
  #画连线
  for(i in 1:nrow(x)){
    circos.link(sector.index1 = x$Gene_1[i], 
                point1 = c(x[i, 4], x[i, 5]),
                sector.index2 = x$Gene_2[i], 
                point2 = c(x[i, 6], x[i, 7]),
                col = paste(x$color[i], "C9", sep = ""), 
                border = FALSE, 
                rou = 0.7
    )}
  
  #画图例
  i <- seq(0,0.995,0.005)
  rect(-1+i/2, #xleft
       -1, #ybottom
       -0.9975+i/2, #xright
       -0.96, #ytop
       col = paste(as.character(color[,1]), "FF", sep = ""),
       border = paste(as.character(color[,1]), "FF", sep = ""))
  text(-0.97, -1.03, "-1")
  text(-0.51, -1.03, "1")
  
}



for (d in 1:32) {
  inputtemp<-immunedata %>%
    filter(Cancer == unique(immunedata$Cancer)[d]  ) %>%
    dplyr::select(gene,all_of(value)) %>%
    na.omit()
  #计算相关性系数
  genecorl <- lapply(colnames(inputtemp),function(x){
    ddd <- genecor.parallel(data = t(inputtemp), cl=1, gene=x) #一定要注意cl参数根据自己电脑cpu线程调整
    ddd  
  })
  genecor <- do.call(rbind, genecorl)
  genecor$cor[is.na(genecor$cor)] <- 0
  genecor$p.value[is.na(genecor$p.value)] <- 1
  write.csv(genecor,file = paste( gene, as.character(unique(immunedata$Cancer)[d]),'免疫相关性分析结果.csv',sep = ' '),row.names =F)
  # 删掉p value = 0的，也就是自己跟自己配对,和相关系数小于0.15的
  
  genecorr1<- genecor %>%
    filter(Gene_1 == gene  ) %>%
    filter(cor > 0.15 | cor < -0.15)
  if(nrow(genecorr1) > 2 ){
    #选取这些内容从新做相关性分析
    inputtemp<-inputtemp[,genecorr1$Gene_2]
    #计算相关性系数
    genecorl <- lapply(colnames(inputtemp),function(x){
      ddd <- genecor.parallel(data = t(inputtemp), cl=1, gene=x) #一定要注意cl参数根据自己电脑cpu线程调整
      ddd  
    })
    genecor <- do.call(rbind, genecorl)
    # 删掉p value = 0的，也就是自己跟自己配对
    genecorr <- genecor[-which(genecor$p.value==0),]
    # 删掉A vs B 和 B vs A其中一个
    genecorrr<-genecorr[!duplicated(genecorr$cor),]
    # 保存到文件
    genecorrr$p.value <- NULL
    pdf(paste(gene,as.character(unique(immunedata$Cancer)[d]),'相关性圈图.pdf',sep=' '), width = 5, height = 5)
    genecor_circleplot(genecorrr)
    dev.off()
  }else {
    genecorr <- genecor[-which(genecor$p.value==0),]
    # 删掉A vs B 和 B vs A其中一个
    genecorrr<-genecorr[!duplicated(genecorr$cor),]
    # 保存到文件
    genecorrr$p.value <- NULL
    pdf(paste(gene,as.character(unique(immunedata$Cancer)[d]),'相关性圈图.pdf',sep=' '), width = 5, height = 5)
    genecor_circleplot1(genecorrr)
    dev.off()
  }
  #免疫浸润相关性线图
  inputtemp<-immunedata %>%
    filter(Cancer == unique(immunedata$Cancer)[d]  ) %>%
    dplyr::select(gene,all_of(value)) %>%
    na.omit()
  for (b in 1:26) {
    test <- cor.test(inputtemp[,gene],inputtemp[,value[b]],exact=FALSE)
    paste(paste0("n = ",length(inputtemp[,gene])),
          paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
          paste0("p.value= ",round(test$p.value,4)),
          sep = ", ")
    p<-ggplot(inputtemp,aes(get(gene),get(value[b])))+
      geom_point(col="#984ea3")+
      geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
      geom_rug(col="#7fc97f")+
      theme_minimal()+
      xlab(paste(gene, 'expression log2(TPM+0.001)' ,sep = ' '    ))+
      ylab(paste(value[b], 'Infiltration Score' ,sep = ' '    ))+
      ## 依靠函数来生成title
      labs(title = paste(as.character(unique(immunedata$Cancer)[d]),
                         paste0("n = ",length(inputtemp[,gene])),
                         paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
                         paste0("p.value= ",round(test$p.value,4)),
                         sep = ", "))+
      theme(plot.title = element_text(hjust = 0.5),
            plot.margin = margin(1, 1, 1, 1, "cm"))
    ggsave(p,filename = paste(gene,value[b],as.character(unique(immunedata$Cancer)[d]),'相关性点图.pdf',sep=' '),width = 4.94,height = 4.72)
  }
  inputtemp<-immunedata %>%
    filter(Cancer == unique(immunedata$Cancer)[d]  ) %>%
    dplyr::select(gene,all_of(value)) %>%
    na.omit()
  med.exp<-median(inputtemp[,gene])
  more.med.exp.index<-which(inputtemp[,gene]>=med.exp)
  less.med.exp.index<-which(inputtemp[,gene]<med.exp)
  inputtemp[,gene][more.med.exp.index]<-paste0('High expression')
  inputtemp[,gene][less.med.exp.index]<-paste0('Low expression')
  inputtemp1<- inputtemp %>% 
    gather("cell_type", "value",-gene)
  inputtemp1[,gene]<-factor(inputtemp1[,gene],levels = c('High expression','Low expression'))
  for (e in 1:26) {
    aaa<-dplyr::filter(inputtemp1,cell_type == unique(inputtemp1$cell_type)[[e]])
    p<-ggboxplot(aaa, x = gene , y = 'value',
                 fill = gene, legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T,scales = "free_y")+
      theme(legend.position='none')+
      ylab(label = 'Infiltration Score')+
      xlab(label = gene )+
      stat_compare_means(method = 't.test',aes(label = ..p.signif..),label.x = 1.5)+
      ggtitle(unique(inputtemp1$cell_type)[[e]])+theme(plot.title = element_text(hjust = 0.5))  #添加标题并居中
    ggsave(p,filename = paste(gene,as.character(unique(immunedata$Cancer)[d]),  unique(inputtemp1$cell_type)[[e]],   '差异表达.pdf',sep=' '),width = 8,height = 6)
  }
  p<-ggboxplot(inputtemp1, x = gene , y = 'value',
               fill = gene, legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T, facet.by = "cell_type", scales = "free_y")+
    theme(legend.position='none')+
    ylab(label = 'Infiltration Score')+
    xlab(label = gene )+
    stat_compare_means(method = 't.test',aes(label = ..p.signif..),label.x = 1.5)
  ggsave(p,filename = paste('合并图', gene,as.character(unique(immunedata$Cancer)[d]), '差异表达.pdf',sep=' '),width = 25,height = 16)
}

Index<-c( 'ACC', 'BLCA' ,'BRCA' ,'CESC', 'CHOL','COAD','READ','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LGG',
          'LIHC','LUAD', 'LUSC','MESO','OV','PAAD','PCPG', 'PRAD','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM'  )
cor_data_df1<-data.frame(symbol= character(),correlation= numeric(),pvalue= numeric(),type= character())
for (i in 1:32) {
  aaa<-immunedata%>%
    filter(Cancer == Index[[i]], Type == 'Tumor')%>%
    dplyr::select(all_of(gene),value)
  y <- as.numeric(aaa[,gene])#开始相关性分析
  cor_data_df <- data.frame(colnames(aaa))
  for (a in 1:27){
    test <- cor.test(as.numeric(aaa[,a]),y,method = "pearson")
    cor_data_df[a,2] <- test$estimate
    cor_data_df[a,3] <- test$p.value
    cor_data_df[a,4] <- Index[[i]]
  }
  names(cor_data_df) <- c("symbol","correlation","pvalue","type")
  cor_data_df$correlation[which(is.na(cor_data_df$correlation))]<-0
  cor_data_df$pvalue[which(is.na(cor_data_df$pvalue))]<-1
  cor_data_df1<-rbind(cor_data_df1,cor_data_df)
}
write.csv (cor_data_df1,'基因与免疫细胞相关性分析结果.csv' , row.names =FALSE)#将文件导出
aaaaa<-colnames(aaa)[2:27]
cor_data_df1<-dplyr::filter(cor_data_df1,symbol != gene)
cor_data_df1$label<-NA
cor_data_df1[cor_data_df1[,"pvalue"]<=0.0001,"label"]<-"****"
cor_data_df1[cor_data_df1[,"pvalue"]<= 0.001 & cor_data_df1[,"pvalue"]>0.0001,"label"]<-"***"
cor_data_df1[cor_data_df1[,"pvalue"]<= 0.01 & cor_data_df1[,"pvalue"]>0.001,"label"]<-"**"
cor_data_df1[cor_data_df1[,"pvalue"]<= 0.05 & cor_data_df1[,"pvalue"]>0.01,"label"]<-"*"
meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:26) {
  meandata1<-dplyr::filter(cor_data_df1,symbol == aaaaa[[i]] )
  meandata[i,1]<-aaaaa[[i]]
  meandata[i,2]<-mean(meandata1$correlation)
}
bbbbb<-as.character(arrange(meandata,meanvalue)$symbol)
cor_data_df1$symbol<-factor(cor_data_df1$symbol,levels = bbbbb)
p<-ggplot(cor_data_df1,aes(x=type,y=symbol))+geom_tile(aes(fill=correlation))+
  scale_fill_gradientn(colors=c("blue","white","red"),guide="colorbar")+
  theme_void()+
  theme(axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))+
  theme(axis.text.y=element_text(angle=0,hjust=1))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_text(aes(label = label),col='black',cex=2.5)
ggsave(p,filename = '基因与免疫细胞相关性分析热图.pdf',width =10,height =8)
#----------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------
#基因与免疫相关基因相关性分析
Index<-c( 'ACC', 'BLCA' ,'BRCA' ,'CESC', 'CHOL','COAD','READ','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG',
          'LIHC','LUAD', 'LUSC','MESO','OV','PAAD','PCPG', 'PRAD','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM'  )
immunegene<-read.table('immunegene.txt',header = T,sep = '\t')
othergene<-read.table('新的基因集.txt',header = T,sep = '\t')
# 计算相关系数的函数
genecor.parallel <- function(data,gene,cl){
  cl <- makeCluster(cl)
  y <- as.numeric(data[gene,])
  rownames <- rownames(data)
  dataframes <- do.call(rbind, parLapply(cl=cl,rownames, function(x){
    dd  <- cor.test(as.numeric(data[x,]), y, type="spearman")
    data.frame(Gene_1=gene, Gene_2=x, cor=dd$estimate, p.value=dd$p.value)
  }))
  stopCluster(cl)
  return(dataframes)
}
genecor_circleplot <- function(x){
  Corr <- data.frame(rbind(data.frame(Gene=x[,1], Correlation=x[,3]), 
                           data.frame(Gene=x[,2], Correlation=x[,3])), stringsAsFactors = F)      
  Corr$Index <- seq(1,nrow(Corr),1) #记录基因的原始排序，记录到Index列
  Corr <- Corr[order(Corr[,1]),] #按照基因名排序
  corrsp <- split(Corr,Corr$Gene)
  corrspe <- lapply(corrsp, function(x){x$Gene_Start<-0
  
  #依次计算每个基因的相关系数总和，作为基因终止位点
  if (nrow(x)==1){x$Gene_End<-1}else{
    x$Gene_End<-sum(abs(x$Correlation))} 
  x})
  GeneID <- do.call(rbind,corrspe)
  GeneID <- GeneID[!duplicated(GeneID$Gene),]
  
  #基因配色
  mycol <- pal_d3("category20c")(20)
  n <- nrow(GeneID)
  GeneID$Color <- mycol[1:n]
  
  #连线的宽度是相关系数的绝对值
  Corr[,2] <- abs(Corr[,2]) 
  corrsl <- split(Corr,Corr$Gene)
  aaaaa <- c()
  corrspl <- lapply(corrsl,function(x){nn<-nrow(x)
  for (i in 1:nn){
    aaaaa[1] <- 0
    aaaaa[i+1] <- x$Correlation[i]+aaaaa[i]}
  bbbbb <- data.frame(V4=aaaaa[1:nn],V5=aaaaa[2:(nn+1)])
  bbbbbb <- cbind(x,bbbbb)
  bbbbbb
  })
  Corr <- do.call(rbind,corrspl)
  
  #根据Index列，把基因恢复到原始排序
  Corr <- Corr[order(Corr$Index),]
  
  #V4是起始位置，V5是终止位置
  #把它写入Links里，start_1和end_1对应Gene_1，start_2和end_2对应Gene_2
  x$start_1 <- Corr$V4[1:(nrow(Corr)/2)]
  x$end_1 <- Corr$V5[1:(nrow(Corr)/2)]
  x$start_2 <- Corr$V4[(nrow(Corr)/2 + 1):nrow(Corr)]
  x$end_2 <- Corr$V5[(nrow(Corr)/2 + 1):nrow(Corr)]
  
  #连线（相关系数）的配色
  #相关系数最大为1，最小-1，此处设置201个颜色
  #-1到0就是前100，0到1就是后100
  color <- data.frame(colorRampPalette(c("#67BE54", "#FFFFFF", "#F82C2B"))(201))
  #根据相关系数的数值，给出相应的颜色
  for (i in 1:nrow(x)){
    x[i,8] <- substring(color[x[i,3] * 100 + 101, 1], 1, 7)
  }
  names(x)[8] <- "color"
  
  #绘图区设置
  #par(mar=rep(0,4))
  circos.clear()
  circos.par(start.degree = 90, #从哪里开始画，沿着逆时针顺序
             gap.degree = 5, #基因bar之间的间隔大小
             track.margin = c(0,0.23), #值越大，基因跟连线的间隔越小
             cell.padding = c(0,0,0,0)
  )
  circos.initialize(factors = GeneID$Gene,
                    xlim = cbind(GeneID$Gene_Start, GeneID$Gene_End))
  
  #先画基因
  circos.trackPlotRegion(ylim = c(0, 1), factors = GeneID$Gene, 
                         track.height = 0.05, #基因线条的胖瘦
                         panel.fun = function(x, y) {
                           name = get.cell.meta.data("sector.index") 
                           i = get.cell.meta.data("sector.numeric.index") 
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           circos.text(x = mean(xlim), y = 1,
                                       labels = name,
                                       cex = 1, #基因ID文字大小
                                       niceFacing = TRUE, #保持基因名的头朝上
                                       facing = "bending", #基因名沿着圆弧方向，还可以是reverse.clockwise
                                       adj = c(0.5, -2.8), #基因名所在位置，分别控制左右和上下
                                       font = 2 #加粗
                           )
                           circos.rect(xleft = xlim[1], 
                                       ybottom = ylim[1],
                                       xright = xlim[2], 
                                       ytop = ylim[2],
                                       col = GeneID$Color[i],
                                       border = GeneID$Color[i])
                           
                           circos.axis(labels.cex = 0.7, 
                                       direction = "outside"
                           )})
  
  #画连线
  for(i in 1:nrow(x)){
    circos.link(sector.index1 = x$Gene_1[i], 
                point1 = c(x[i, 4], x[i, 5]),
                sector.index2 = x$Gene_2[i], 
                point2 = c(x[i, 6], x[i, 7]),
                col = paste(x$color[i], "C9", sep = ""), 
                border = FALSE, 
                rou = 0.7
    )}
  
  #画图例
  i <- seq(0,0.995,0.005)
  rect(-1+i/2, #xleft
       -1, #ybottom
       -0.9975+i/2, #xright
       -0.96, #ytop
       col = paste(as.character(color[,1]), "FF", sep = ""),
       border = paste(as.character(color[,1]), "FF", sep = ""))
  text(-0.97, -1.03, "-1")
  text(-0.51, -1.03, "1")
}
cor_data_df1<-data.frame(symbol= character(),correlation= numeric(),pvalue= numeric(),type= character())
if(gene %in% immunegene$symbol){
  for (i in 1:33) {
    aaa<-drawdata%>%
      filter(Cancer == Index[[i]], Type == 'Tumor')%>%
      dplyr::select(all_of(gene),immunegene$symbol)
    y <- as.numeric(aaa[,gene])#开始相关性分析
    cor_data_df <- data.frame(colnames(aaa))
    for (a in 1:221){
      test <- cor.test(as.numeric(aaa[,a]),y,method = "pearson")
      cor_data_df[a,2] <- test$estimate
      cor_data_df[a,3] <- test$p.value
      cor_data_df[a,4] <- Index[[i]]
    }
    names(cor_data_df) <- c("symbol","correlation","pvalue","type")
    cor_data_df1<-rbind(cor_data_df1,cor_data_df)
  }
}else {
  for (i in 1:33) {
    aaa<-drawdata%>%
      filter(Cancer == Index[[i]], Type == 'Tumor')%>%
      dplyr::select(all_of(gene),immunegene$symbol)
    y <- as.numeric(aaa[,gene])#开始相关性分析
    cor_data_df <- data.frame(colnames(aaa))
    for (a in 1:222){
      test <- cor.test(as.numeric(aaa[,a]),y,method = "pearson")
      cor_data_df[a,2] <- test$estimate
      cor_data_df[a,3] <- test$p.value
      cor_data_df[a,4] <- Index[[i]]
    }
    names(cor_data_df) <- c("symbol","correlation","pvalue","type")
    cor_data_df1<-rbind(cor_data_df1,cor_data_df)
  }
}

write.csv (cor_data_df1,'基因与免疫相关基因相关性分析结果.csv' , row.names =FALSE)#将文件导出

xx<-filter(immunegene,type == 'chemokine')$symbol
chemokine<-cor_data_df1%>%
  filter(symbol  %in%   xx )
chemokine<-na.omit(chemokine)#比如某个基因在某个肿瘤中表达值同位一个数值，会导致相关性结果为NA
chemokine$label<-NA
chemokine[chemokine[,"pvalue"]<=0.0001,"label"]<-"****"
chemokine[chemokine[,"pvalue"]<= 0.001 & chemokine[,"pvalue"]>0.0001,"label"]<-"***"
chemokine[chemokine[,"pvalue"]<= 0.01 & chemokine[,"pvalue"]>0.001,"label"]<-"**"
chemokine[chemokine[,"pvalue"]<= 0.05 & chemokine[,"pvalue"]>0.01,"label"]<-"*"

aaaaa<-unique(chemokine$symbol)
meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(aaaaa)) {
  meandata1<-dplyr::filter(chemokine,symbol == aaaaa[[i]] )
  meandata[i,1]<-aaaaa[[i]]
  meandata[i,2]<-mean(meandata1$correlation)
}
bbbbb<-as.character(arrange(meandata,meanvalue)$symbol)
chemokine$symbol<-factor(chemokine$symbol,levels = bbbbb)

meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(unique(chemokine$type))) {
  meandata2<-dplyr::filter(chemokine,type == unique(chemokine$type)[[i]] )
  meandata[i,1]<-Index[[i]]
  meandata[i,2]<-mean(meandata2$correlation)
}
ccccc<-as.character(arrange(meandata,-meanvalue)$symbol)
chemokine$type<-factor(chemokine$type,levels = ccccc)


xx<-filter(immunegene,type == 'receptor')$symbol
receptor<-cor_data_df1%>%
  filter(symbol  %in%   xx )
receptor<-na.omit(receptor)#比如某个基因在某个肿瘤中表达值同位一个数值，会导致相关性结果为NA
receptor$label<-NA
receptor[receptor[,"pvalue"]<=0.0001,"label"]<-"****"
receptor[receptor[,"pvalue"]<= 0.001 & receptor[,"pvalue"]>0.0001,"label"]<-"***"
receptor[receptor[,"pvalue"]<= 0.01 & receptor[,"pvalue"]>0.001,"label"]<-"**"
receptor[receptor[,"pvalue"]<= 0.05 & receptor[,"pvalue"]>0.01,"label"]<-"*"

aaaaa<-unique(receptor$symbol)
meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(aaaaa)) {
  meandata1<-dplyr::filter(receptor,symbol == aaaaa[[i]] )
  meandata[i,1]<-aaaaa[[i]]
  meandata[i,2]<-mean(meandata1$correlation)
}
bbbbb<-as.character(arrange(meandata,meanvalue)$symbol)
receptor$symbol<-factor(receptor$symbol,levels = bbbbb)

meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(unique(receptor$type))) {
  meandata2<-dplyr::filter(receptor,type == unique(receptor$type)[[i]] )
  meandata[i,1]<-Index[[i]]
  meandata[i,2]<-mean(meandata2$correlation)
}
ccccc<-as.character(arrange(meandata,-meanvalue)$symbol)
receptor$type<-factor(receptor$type,levels = ccccc)



xx<-filter(immunegene,type == 'MHC')$symbol
MHC<-cor_data_df1%>%
  filter(symbol  %in%   xx )
MHC<-na.omit(MHC)#比如某个基因在某个肿瘤中表达值同位一个数值，会导致相关性结果为NA
MHC$label<-NA
MHC[MHC[,"pvalue"]<=0.0001,"label"]<-"****"
MHC[MHC[,"pvalue"]<= 0.001 & MHC[,"pvalue"]>0.0001,"label"]<-"***"
MHC[MHC[,"pvalue"]<= 0.01 & MHC[,"pvalue"]>0.001,"label"]<-"**"
MHC[MHC[,"pvalue"]<= 0.05 & MHC[,"pvalue"]>0.01,"label"]<-"*"

aaaaa<-unique(MHC$symbol)
meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(aaaaa)) {
  meandata1<-dplyr::filter(MHC,symbol == aaaaa[[i]] )
  meandata[i,1]<-aaaaa[[i]]
  meandata[i,2]<-mean(meandata1$correlation)
}
bbbbb<-as.character(arrange(meandata,meanvalue)$symbol)
MHC$symbol<-factor(MHC$symbol,levels = bbbbb)

meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(unique(MHC$type))) {
  meandata2<-dplyr::filter(MHC,type == unique(MHC$type)[[i]] )
  meandata[i,1]<-Index[[i]]
  meandata[i,2]<-mean(meandata2$correlation)
}
ccccc<-as.character(arrange(meandata,-meanvalue)$symbol)
MHC$type<-factor(MHC$type,levels = ccccc)



xx<-filter(immunegene,type == 'Immunoinhibitor')$symbol
Immunoinhibitor<-cor_data_df1%>%
  filter(symbol  %in%   xx )
Immunoinhibitor<-na.omit(Immunoinhibitor)#比如某个基因在某个肿瘤中表达值同位一个数值，会导致相关性结果为NA
Immunoinhibitor$label<-NA
Immunoinhibitor[Immunoinhibitor[,"pvalue"]<=0.0001,"label"]<-"****"
Immunoinhibitor[Immunoinhibitor[,"pvalue"]<= 0.001 & Immunoinhibitor[,"pvalue"]>0.0001,"label"]<-"***"
Immunoinhibitor[Immunoinhibitor[,"pvalue"]<= 0.01 & Immunoinhibitor[,"pvalue"]>0.001,"label"]<-"**"
Immunoinhibitor[Immunoinhibitor[,"pvalue"]<= 0.05 & Immunoinhibitor[,"pvalue"]>0.01,"label"]<-"*"

aaaaa<-unique(Immunoinhibitor$symbol)
meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(aaaaa)) {
  meandata1<-dplyr::filter(Immunoinhibitor,symbol == aaaaa[[i]] )
  meandata[i,1]<-aaaaa[[i]]
  meandata[i,2]<-mean(meandata1$correlation)
}
bbbbb<-as.character(arrange(meandata,meanvalue)$symbol)
Immunoinhibitor$symbol<-factor(Immunoinhibitor$symbol,levels = bbbbb)

meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(unique(Immunoinhibitor$type))) {
  meandata2<-dplyr::filter(Immunoinhibitor,type == unique(Immunoinhibitor$type)[[i]] )
  meandata[i,1]<-Index[[i]]
  meandata[i,2]<-mean(meandata2$correlation)
}
ccccc<-as.character(arrange(meandata,-meanvalue)$symbol)
Immunoinhibitor$type<-factor(Immunoinhibitor$type,levels = ccccc)



xx<-filter(immunegene,type == 'Immunostimulator')$symbol
Immunostimulator<-cor_data_df1%>%
  filter(symbol  %in%   xx )
Immunostimulator<-na.omit(Immunostimulator)#比如某个基因在某个肿瘤中表达值同位一个数值，会导致相关性结果为NA
Immunostimulator$label<-NA
Immunostimulator[Immunostimulator[,"pvalue"]<=0.0001,"label"]<-"****"
Immunostimulator[Immunostimulator[,"pvalue"]<= 0.001 & Immunostimulator[,"pvalue"]>0.0001,"label"]<-"***"
Immunostimulator[Immunostimulator[,"pvalue"]<= 0.01 & Immunostimulator[,"pvalue"]>0.001,"label"]<-"**"
Immunostimulator[Immunostimulator[,"pvalue"]<= 0.05 & Immunostimulator[,"pvalue"]>0.01,"label"]<-"*"

aaaaa<-unique(Immunostimulator$symbol)
meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(aaaaa)) {
  meandata1<-dplyr::filter(Immunostimulator,symbol == aaaaa[[i]] )
  meandata[i,1]<-aaaaa[[i]]
  meandata[i,2]<-mean(meandata1$correlation)
}
bbbbb<-as.character(arrange(meandata,meanvalue)$symbol)
Immunostimulator$symbol<-factor(Immunostimulator$symbol,levels = bbbbb)

meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(unique(Immunostimulator$type))) {
  meandata2<-dplyr::filter(Immunostimulator,type == unique(Immunostimulator$type)[[i]] )
  meandata[i,1]<-Index[[i]]
  meandata[i,2]<-mean(meandata2$correlation)
}
ccccc<-as.character(arrange(meandata,-meanvalue)$symbol)
Immunostimulator$type<-factor(Immunostimulator$type,levels = ccccc)



xx<-filter(immunegene,type == 'm6a')$symbol
m6a<-cor_data_df1%>%
  filter(symbol  %in%   xx )
m6a<-na.omit(m6a)#比如某个基因在某个肿瘤中表达值同位一个数值，会导致相关性结果为NA
m6a$label<-NA
m6a[m6a[,"pvalue"]<=0.0001,"label"]<-"****"
m6a[m6a[,"pvalue"]<= 0.001 & m6a[,"pvalue"]>0.0001,"label"]<-"***"
m6a[m6a[,"pvalue"]<= 0.01 & m6a[,"pvalue"]>0.001,"label"]<-"**"
m6a[m6a[,"pvalue"]<= 0.05 & m6a[,"pvalue"]>0.01,"label"]<-"*"

aaaaa<-unique(m6a$symbol)
meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(aaaaa)) {
  meandata1<-dplyr::filter(m6a,symbol == aaaaa[[i]] )
  meandata[i,1]<-aaaaa[[i]]
  meandata[i,2]<-mean(meandata1$correlation)
}
bbbbb<-as.character(arrange(meandata,meanvalue)$symbol)
m6a$symbol<-factor(m6a$symbol,levels = bbbbb)

meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(unique(m6a$type))) {
  meandata2<-dplyr::filter(m6a,type == unique(m6a$type)[[i]] )
  meandata[i,1]<-Index[[i]]
  meandata[i,2]<-mean(meandata2$correlation)
}
ccccc<-as.character(arrange(meandata,-meanvalue)$symbol)
m6a$type<-factor(m6a$type,levels = ccccc)



xx<-filter(immunegene,type == 'Ferrotosis')$symbol
Ferrotosis<-cor_data_df1%>%
  filter(symbol  %in%   xx )
Ferrotosis<-na.omit(Ferrotosis)#比如某个基因在某个肿瘤中表达值同位一个数值，会导致相关性结果为NA
Ferrotosis$label<-NA
Ferrotosis[Ferrotosis[,"pvalue"]<=0.0001,"label"]<-"****"
Ferrotosis[Ferrotosis[,"pvalue"]<= 0.001 & Ferrotosis[,"pvalue"]>0.0001,"label"]<-"***"
Ferrotosis[Ferrotosis[,"pvalue"]<= 0.01 & Ferrotosis[,"pvalue"]>0.001,"label"]<-"**"
Ferrotosis[Ferrotosis[,"pvalue"]<= 0.05 & Ferrotosis[,"pvalue"]>0.01,"label"]<-"*"

aaaaa<-unique(Ferrotosis$symbol)
meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(aaaaa)) {
  meandata1<-dplyr::filter(Ferrotosis,symbol == aaaaa[[i]] )
  meandata[i,1]<-aaaaa[[i]]
  meandata[i,2]<-mean(meandata1$correlation)
}
bbbbb<-as.character(arrange(meandata,meanvalue)$symbol)
Ferrotosis$symbol<-factor(Ferrotosis$symbol,levels = bbbbb)

meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(unique(Ferrotosis$type))) {
  meandata2<-dplyr::filter(Ferrotosis,type == unique(Ferrotosis$type)[[i]] )
  meandata[i,1]<-Index[[i]]
  meandata[i,2]<-mean(meandata2$correlation)
}
ccccc<-as.character(arrange(meandata,-meanvalue)$symbol)
Ferrotosis$type<-factor(Ferrotosis$type,levels = ccccc)



p<-ggplot(m6a,aes(x=type,y=symbol))+geom_tile(aes(fill=correlation))+
  scale_fill_gradientn(colors=c("blue","white","red"),guide="colorbar")+
  theme_void()+
  theme(axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))+
  theme(axis.text.y=element_text(angle=0,hjust=1))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_text(aes(label = label),col='black',cex=2.5)
ggsave(p,filename = '基因与m6a相关性分析热图.pdf',width =10,height =4)

p<-ggplot(Ferrotosis,aes(x=type,y=symbol))+geom_tile(aes(fill=correlation))+
  scale_fill_gradientn(colors=c("blue","white","red"),guide="colorbar")+
  theme_void()+
  theme(axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))+
  theme(axis.text.y=element_text(angle=0,hjust=1))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_text(aes(label = label),col='black',cex=2.5)
ggsave(p,filename = '基因与铁死亡相关性分析热图.pdf',width =10,height =20)

p<-ggplot(chemokine,aes(x=type,y=symbol))+geom_tile(aes(fill=correlation))+
  scale_fill_gradientn(colors=c("blue","white","red"),guide="colorbar")+
  theme_void()+
  theme(axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))+
  theme(axis.text.y=element_text(angle=0,hjust=1))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_text(aes(label = label),col='black',cex=2.5)
ggsave(p,filename = '基因与趋化因子相关性分析热图.pdf',width =8,height =10)

p<-ggplot(receptor,aes(x=type,y=symbol))+geom_tile(aes(fill=correlation))+
  scale_fill_gradientn(colors=c("blue","white","red"),guide="colorbar")+
  theme_void()+
  theme(axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))+
  theme(axis.text.y=element_text(angle=0,hjust=1))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_text(aes(label = label),col='black',cex=2.5)
ggsave(p,filename = '基因与趋化因子受体相关性分析热图.pdf',width =8,height = 4.72)

p<-ggplot(MHC,aes(x=type,y=symbol))+geom_tile(aes(fill=correlation))+
  scale_fill_gradientn(colors=c("blue","white","red"),guide="colorbar")+
  theme_void()+
  theme(axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))+
  theme(axis.text.y=element_text(angle=0,hjust=1))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_text(aes(label = label),col='black',cex=2.5)
ggsave(p,filename = '基因与MHC基因相关性分析热图.pdf',width =8,height = 4.72)

p<-ggplot(Immunoinhibitor,aes(x=type,y=symbol))+geom_tile(aes(fill=correlation))+
  scale_fill_gradientn(colors=c("blue","white","red"),guide="colorbar")+
  theme_void()+
  theme(axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))+
  theme(axis.text.y=element_text(angle=0,hjust=1))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_text(aes(label = label),col='black',cex=2.5)
ggsave(p,filename = '基因与免疫抑制基因相关性分析热图.pdf',width =8,height = 6)

p<-ggplot(Immunostimulator,aes(x=type,y=symbol))+geom_tile(aes(fill=correlation))+
  scale_fill_gradientn(colors=c("blue","white","red"),guide="colorbar")+
  theme_void()+
  theme(axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))+
  theme(axis.text.y=element_text(angle=0,hjust=1))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_text(aes(label = label),col='black',cex=2.5)
ggsave(p,filename = '基因与免疫激活基因相关性分析热图.pdf',width =8,height = 10)


####单因素森林图
genename<-gene
#OS
Index<-sort(unique(OSdata$Cancer))
outTab=data.frame()
for (i in 1:length(Index)) {
  dataa<-filter(OSdata,Cancer ==Index[i] )
  cox <- coxph(Surv(time, status) ~ dataa[,genename], data = dataa)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=Index[i],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="OS单因素结果.xls",sep="\t",row.names=F,quote=F)
rt=read.table("OS单因素结果.xls",header=T,sep="\t",row.names=1,check.names=F)
#把含有0或者NA的行全都去掉
rt[rt == 0] <- NA
rt<-na.omit(rt)
rt$HR<-log2(rt$HR)
rt$HR.95L<-log2(rt$HR.95L)
rt$HR.95H<-log2(rt$HR.95H)
rt<-arrange(rt,pvalue)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue") 
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   
pdf(file="OS_forest.pdf",
    width = 6,            
    height = 8,           
)
forestplot(tabletext, 
           zero = 0,
           lwd.zero = 2,
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=F,
           lwd.ci=2,
           boxsize=0.4,
           xlab="log2(Hazard ratio)"
)
dev.off()

#DFI
DFI_DATA1<-filter(DFI_DATA,Cancer != "GBM")
Index<-sort(unique(DFI_DATA1$Cancer))
outTab=data.frame()
for (i in 1:length(Index)) {
  dataa<-filter(DFI_DATA1,Cancer ==Index[i] )
  cox <- coxph(Surv(DFI.time, DFI) ~ dataa[,genename], data = dataa)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=Index[i],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="DFI单因素结果.xls",sep="\t",row.names=F,quote=F)
rt=read.table("DFI单因素结果.xls",header=T,sep="\t",row.names=1,check.names=F)
#把含有0或者NA的行全都去掉
rt[rt == 0] <- NA
rt<-na.omit(rt)
rt$HR<-log2(rt$HR)
rt$HR.95L<-log2(rt$HR.95L)
rt$HR.95H<-log2(rt$HR.95H)
rt<-arrange(rt,pvalue)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue") 
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   
pdf(file="DFI_forest.pdf",
    width = 6,            
    height = 8,           
)
forestplot(tabletext, 
           zero = 0,
           lwd.zero = 2,
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=F,
           lwd.ci=2,
           boxsize=0.4,
           xlab="log2(Hazard ratio)"
)
dev.off()
#DSS
Index<-sort(unique(DSS_DATA$Cancer))
outTab=data.frame()
for (i in 1:length(Index)) {
  dataa<-filter(DSS_DATA,Cancer ==Index[i] )
  cox <- coxph(Surv(DSS.time, DSS) ~ dataa[,genename], data = dataa)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=Index[i],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="DSS单因素结果.xls",sep="\t",row.names=F,quote=F)
rt=read.table("DSS单因素结果.xls",header=T,sep="\t",row.names=1,check.names=F)
#把含有0或者NA的行全都去掉
rt[rt == 0] <- NA
rt<-na.omit(rt)
rt$HR<-log2(rt$HR)
rt$HR.95L<-log2(rt$HR.95L)
rt$HR.95H<-log2(rt$HR.95H)
rt<-arrange(rt,pvalue)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue") 
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   
pdf(file="DSS_forest.pdf",
    width = 6,            
    height = 8,           
)
forestplot(tabletext, 
           zero = 0,
           lwd.zero = 2,
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=F,
           lwd.ci=2,
           boxsize=0.4,
           xlab="log2(Hazard ratio)"
)
dev.off()
###PFI
Index<-sort(unique(PFI_DATA$Cancer))
outTab=data.frame()
for (i in 1:length(Index)) {
  dataa<-filter(PFI_DATA,Cancer ==Index[i] )
  cox <- coxph(Surv(PFI.time, PFI) ~ dataa[,genename], data = dataa)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=Index[i],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="PFI单因素结果.xls",sep="\t",row.names=F,quote=F)
rt=read.table("PFI单因素结果.xls",header=T,sep="\t",row.names=1,check.names=F)
#把含有0或者NA的行全都去掉
rt[rt == 0] <- NA
rt<-na.omit(rt)
rt$HR<-log2(rt$HR)
rt$HR.95L<-log2(rt$HR.95L)
rt$HR.95H<-log2(rt$HR.95H)
rt<-arrange(rt,pvalue)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue") 
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   
pdf(file="PFI_forest.pdf",
    width = 6,            
    height = 8,           
)
forestplot(tabletext, 
           zero = 0,
           lwd.zero = 2,
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=F,
           lwd.ci=2,
           boxsize=0.4,
           xlab="log2(Hazard ratio)"
)
dev.off()
#------------------------------------------------------------------



#----------------------------------------------------------
#TME分析1
Index<-c( 'ACC', 'BLCA' ,'BRCA' ,'CESC', 'CHOL','COAD','READ','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG',
          'LIHC','LUAD', 'LUSC','MESO','OV','PAAD','PCPG', 'PRAD','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM'  )
TMEscores1<-TMEscores
TMEscores1<-TMEscores1[,-16]
TMEscores1<-TMEscores1[,-6]
TMEscores1<-TMEscores1[,-2]
for (i in 1:33) {
  pdata<-TMEscores1 %>%
    dplyr::filter(Cancer == Index[[i]])%>%
    dplyr::select(colnames(TMEscores1)[1:13],all_of(gene))
  pdata$group<-ifelse(pdata[,gene]>= median(pdata[,gene]),"High expression","Low expression")
  pdata<-pdata[,-which(colnames(pdata) == gene)]
  pdata_melt <- reshape2::melt(pdata,
                               id.vars = c("ID","group"),
                               variable.name ="Signature",
                               value.name = "Signature_score")
  #开始画图
  # 使用ggplo2画图
  c <- ggplot(pdata_melt,
              aes(x=Signature, y=Signature_score, 
                  fill = group, #按类填充颜色
                  color = group)) + #按类给边框着色
    geom_boxplot(notch = F, alpha = 0.95, 
                 outlier.shape = 16,
                 outlier.colour = "black", #outlier点用黑色
                 outlier.size = 0.65) +
    #自定义配色
    scale_fill_manual(values= c("#E31A1C","#E7B800","#2E9FDF")) +
    #ggtitle(paste0(Index[[i]],"_signature score")) + 
    labs(title=paste0(Index[[i]]," signature score"))+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")+
    theme(plot.title = element_text(hjust = 0.5)) 
  
  # 标注*
  # `****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1
  p <- c + stat_compare_means(label = "p.signif")
  ggsave(p, filename =paste0(Index[[i]],"_TME高低差异图.pdf"), width = 13, height = 6)
  print(paste0(i,"_TME"))
}

#相关性热图
cor_data_df1<-data.frame(symbol= character(),correlation= numeric(),pvalue= numeric(),type= character())
for (L in 1:33) {
  aaa<-TMEscores1 %>%
    dplyr::filter(Cancer == Index[[L]])%>%
    dplyr::select(colnames(TMEscores1)[1:13],all_of(gene))
  aaa<-aaa[,-1]
  y <- as.numeric(aaa[,gene])#开始相关性分析
  cor_data_df <- data.frame(colnames(aaa))
  for (a in 1:13){
    test <- cor.test(as.numeric(aaa[,a]),y,method = "pearson")
    cor_data_df[a,2] <- test$estimate
    cor_data_df[a,3] <- test$p.value
    cor_data_df[a,4] <- Index[[L]]
  }
  names(cor_data_df) <- c("symbol","correlation","pvalue","type")
  cor_data_df$correlation[which(is.na(cor_data_df$correlation))]<-0
  cor_data_df$pvalue[which(is.na(cor_data_df$pvalue))]<-1
  cor_data_df1<-rbind(cor_data_df1,cor_data_df)
}
write.csv (cor_data_df1,'基因与TME相关性分析结果.csv' , row.names =FALSE)#将文件导出
aaaaa<-colnames(aaa)[1:12]
cor_data_df1<-dplyr::filter(cor_data_df1,symbol != gene)
cor_data_df1$label<-NA
cor_data_df1[cor_data_df1[,"pvalue"]<=0.0001,"label"]<-"****"
cor_data_df1[cor_data_df1[,"pvalue"]<= 0.001 & cor_data_df1[,"pvalue"]>0.0001,"label"]<-"***"
cor_data_df1[cor_data_df1[,"pvalue"]<= 0.01 & cor_data_df1[,"pvalue"]>0.001,"label"]<-"**"
cor_data_df1[cor_data_df1[,"pvalue"]<= 0.05 & cor_data_df1[,"pvalue"]>0.01,"label"]<-"*"
meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:12) {
  meandata1<-dplyr::filter(cor_data_df1,symbol == aaaaa[[i]] )
  meandata[i,1]<-aaaaa[[i]]
  meandata[i,2]<-mean(meandata1$correlation)
}
bbbbb<-as.character(arrange(meandata,meanvalue)$symbol)
cor_data_df1$symbol<-factor(cor_data_df1$symbol,levels = bbbbb)
for (i in 1:33) {
  meandata2<-dplyr::filter(cor_data_df1,type == Index[[i]] )
  meandata[i,1]<-Index[[i]]
  meandata[i,2]<-mean(meandata2$correlation)
}
ccccc<-as.character(arrange(meandata,-meanvalue)$symbol)
cor_data_df1$type<-factor(cor_data_df1$type,levels = ccccc)

p<-ggplot(cor_data_df1,aes(x=type,y=symbol))+geom_tile(aes(fill=correlation))+
  scale_fill_gradientn(colors=c("blue","white","red"),guide="colorbar")+
  theme_void()+
  theme(axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))+
  theme(axis.text.y=element_text(angle=0,hjust=1))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_text(aes(label = label),col='black',cex=2.5)
ggsave(p,filename ="基因与TME相关性分析热图.pdf",width =10,height =6)




#----------------------------------------------------------
#TME分析2
Index<-c( 'ACC', 'BLCA' ,'BRCA' ,'CESC', 'CHOL','COAD','READ','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG',
          'LIHC','LUAD', 'LUSC','MESO','OV','PAAD','PCPG', 'PRAD','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM'  )
Indexx<-c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity")

for (scoretype in Indexx) {
  cor_data_df<-data.frame(correlation= numeric(),pvalue= numeric(),type= character())
  if(gene %in% colnames(estimateScoresdata)){
    for (i in 1:33) {
      aaa<-estimateScoresdata%>%
        filter(Cancer == Index[[i]])
      test <- cor.test(aaa[,gene],aaa[,scoretype],exact=FALSE)
      cor_data_df[i,1] <- test$estimate
      cor_data_df[i,2] <- test$p.value
      cor_data_df[i,3] <- Index[[i]]
      p<-ggplot(aaa,aes(get(gene),get(scoretype)))+
        geom_point(col="#984ea3")+
        geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
        geom_rug(col="#7fc97f")+
        theme_minimal()+
        xlab(paste(gene, 'expression log2(TPM+0.001)' ,sep = ' '    ))+
        ylab(scoretype)+
        ## 依靠函数来生成title
        labs(title = paste(Index[[i]], paste0("n = ",length(aaa[,scoretype])),
                           paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
                           paste0("p.value= ",round(test$p.value,4)),
                           sep = ", "))+
        theme(plot.title = element_text(hjust = 0.5),
              plot.margin = margin(1, 1, 1, 1, "cm"))
      ggsave(p,filename = paste(Index[[i]], gene,"_",scoretype,'相关性图.pdf',sep=' '),width = 4.94,height = 4.72)
    }
    cor_data_df$label<-NA
    cor_data_df[cor_data_df$pvalue<0.05 & cor_data_df$correlation >0,"label"]<-"positive"
    cor_data_df[cor_data_df$pvalue<0.05 & cor_data_df$correlation <0,"label"]<-"negtive"
    cor_data_df[is.na(cor_data_df$label),"label"]<-"non-significant"
    write.csv(cor_data_df,"相关分析结果.csv", row.names =FALSE)
    #棒棒糖
    if(length(unique(cor_data_df$label)) == 3) {
      p<-ggdotchart(cor_data_df, x = "type", y = "correlation",
                    color = "label",                                # Color by groups
                    palette = c("#00AFBB", "grey", "#FC4E07"), # Custom color palette
                    sorting = "descending",                       # Sort value in descending order
                    add = "segments",                             # Add segments from y = 0 to dots
                    add.params = list(color = "lightgray", size = 2), # Change segment color and size
                    #group = "label",                                # Order by groups
                    dot.size = 8,                                 # Large dot size
                    label = round(cor_data_df$correlation,1),                        # Add mpg values as dot labels
                    font.label = list(color = "white", size = 9,  vjust = 0.5),  # Adjust label parameters
                    ggtheme = theme_pubr()  )+
        geom_hline(yintercept = 0, linetype = 2, color = "lightgray")
      ggsave(p,filename = '相关棒棒糖图.pdf',width =12,height =6)
    }else{
      p<-ggdotchart(cor_data_df, x = "type", y = "correlation",
                    color = "label",                                # Color by groups
                    palette = c("#00AFBB",  "#FC4E07"), # Custom color palette
                    sorting = "descending",                       # Sort value in descending order
                    add = "segments",                             # Add segments from y = 0 to dots
                    add.params = list(color = "lightgray", size = 2), # Change segment color and size
                    #group = "label",                                # Order by groups
                    dot.size = 8,                                 # Large dot size
                    label = round(cor_data_df$correlation,1),                        # Add mpg values as dot labels
                    font.label = list(color = "white", size = 9,  vjust = 0.5),  # Adjust label parameters
                    ggtheme = theme_pubr()  )+
        geom_hline(yintercept = 0, linetype = 2, color = "lightgray")
      ggsave(p,filename = '相关棒棒糖图.pdf',width =12,height =6)
    }
    lengthh<-length(cor_data_df$type)
    radardata<-dplyr::select(cor_data_df,type,correlation)
    rownames(radardata)<-radardata$type
    radardata<-as.data.frame(t(radardata))
    radardata<-as.data.frame(radardata[-1,])
    radardata[,1:lengthh]<-apply(radardata[,1:lengthh], 1, as.numeric)
    radardata<-dplyr::select(radardata,colnames(sort(radardata[1,1:lengthh])))
    #radardata[,1:lengthh]<-apply(radardata[,1:lengthh], 1, as.numeric)
    maxValue=ceiling(max(abs(radardata))*10)/10
    radardata=rbind(rep(maxValue,ncol(radardata)),rep(-maxValue,ncol(radardata)),radardata)
    colors="red"
    pdf(file="radar.pdf",height=7,width=7)
    radarchart( radardata, axistype=1 , 
                pcol=colors,                 
                plwd=2 ,                     
                plty=1,
                pty=16,
                cglcol="grey",               
                cglty=1,                     
                caxislabels=seq(-maxValue,maxValue,maxValue/2),
                cglwd=1.2,                  
                axislabcol="blue",          
                vlcex=0.8,
                title=scoretype
    )
    dev.off()
  }
  print(scoretype)
}

cor_data_df1<-data.frame(symbol= character(),correlation= numeric(),pvalue= numeric(),type= character())
for (i in 1:33) {
  aaa<-estimateScoresdata%>%
    filter(Cancer == Index[[i]])%>%
    dplyr::select(all_of(gene),all_of(Indexx))
  y <- as.numeric(aaa[,gene])#开始相关性分析
  cor_data_df <- data.frame(colnames(aaa))
  for (a in 1:5){
    test <- cor.test(as.numeric(aaa[,a]),y,method = "pearson")
    cor_data_df[a,2] <- test$estimate
    cor_data_df[a,3] <- test$p.value
    cor_data_df[a,4] <- Index[[i]]
  }
  names(cor_data_df) <- c("symbol","correlation","pvalue","type")
  cor_data_df1<-rbind(cor_data_df1,cor_data_df)
}

correlationscore<-cor_data_df1%>%
  filter(symbol  %in%   Indexx )
correlationscore<-na.omit(correlationscore)#比如某个基因在某个肿瘤中表达值同位一个数值，会导致相关性结果为NA
correlationscore$label<-NA
correlationscore[correlationscore[,"pvalue"]<=0.0001,"label"]<-"****"
correlationscore[correlationscore[,"pvalue"]<= 0.001 & correlationscore[,"pvalue"]>0.0001,"label"]<-"***"
correlationscore[correlationscore[,"pvalue"]<= 0.01 & correlationscore[,"pvalue"]>0.001,"label"]<-"**"
correlationscore[correlationscore[,"pvalue"]<= 0.05 & correlationscore[,"pvalue"]>0.01,"label"]<-"*"

aaaaa<-unique(correlationscore$symbol)
meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(aaaaa)) {
  meandata1<-dplyr::filter(correlationscore,symbol == aaaaa[[i]] )
  meandata[i,1]<-aaaaa[[i]]
  meandata[i,2]<-mean(meandata1$correlation)
}
bbbbb<-as.character(arrange(meandata,meanvalue)$symbol)
correlationscore$symbol<-factor(correlationscore$symbol,levels = bbbbb)

meandata<-data.frame(symbol= character(),meanvalue= numeric())
for (i in 1:length(unique(correlationscore$type))) {
  meandata2<-dplyr::filter(correlationscore,type == unique(correlationscore$type)[[i]] )
  meandata[i,1]<-Index[[i]]
  meandata[i,2]<-mean(meandata2$correlation)
}
ccccc<-as.character(arrange(meandata,-meanvalue)$symbol)
correlationscore$type<-factor(correlationscore$type,levels = ccccc)

p<-ggplot(correlationscore,aes(x=type,y=symbol))+geom_tile(aes(fill=correlation))+
  scale_fill_gradientn(colors=c("blue","white","red"),guide="colorbar")+
  theme_void()+
  theme(axis.text.x=element_text(angle=90,vjust=0.3,hjust=1))+
  theme(axis.text.y=element_text(angle=0,hjust=1))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_text(aes(label = label),col='black',cex=2.5)
ggsave(p,filename = '相关性分析热图.pdf',width =9,height =2)
