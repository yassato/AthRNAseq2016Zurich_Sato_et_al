########################################################
###Sato et al. Field RNA-Seq of Arabidopsis herbivory###
########################################################

#load data-------------------------------------------------------------
dir.input = "170405_SatoAthRNA"
dir.output = "170405_SatoAthRNA"

fn.description = "160101-1_GeneDescription_Ath"
fn.sample.index.merged = "180104_Sample-IndexPrimer_SatoAth.csv"

rpm.date = "170405"
species = "SatoAth"

fn.rawcnt = sprintf("%s-2_rawcnt_%s", rpm.date, species)
fn.rpm = sprintf("%s-3_rpm_%s", rpm.date, species)

exec.date = "170418_SatoAth"

primerID = read.csv(fn.sample.index.merged, header = T)
load(fn.rawcnt)
load(fn.description)
load(fn.rpm)

load("./180907_GO/ulg.TAIR_180907") #load AGI code with GO terms
#-----------------------------------------------------------------------

#data processing--------------------------------------------------------

#excl. empty samples
primerID = primerID[1:182,]
rpm = rpm[,1:182]

#excl. low rawcnt
hist(log(apply(rawcnt[,1:182],2,sum),2), main = "Lower 5% Limit = 13.2", las=1, breaks=seq(from = 5, to = 25, by=1), xlab="log2(Raw read count)")
rawcnt.l95 = quantile(apply(rawcnt[,1:182],2,sum), 0.05)
abline(v=log(rawcnt.l95,2),lty=2)
excl.id = which(apply(rawcnt[,1:182],2,sum)<rawcnt.l95)

primerID = primerID[-excl.id,]
rpm = rpm[,-excl.id]

#excl. ERCC spikein, virus, and transposon
excl.genes = which(factor(des$Type)!="transposable_element_gene"&factor(des$Type)!="ercc_spikein"&factor(des$Type)!="virus")
rpm = rpm[excl.genes,]
des = des[excl.genes,]

#excl. log(rpm+1)==0
hist(apply(log(rpm+1,2),1,mean),2,breaks=seq(from = 0, to = 16,by=0.25), xlab="Mean_log2(rpm+1)", main="Total 24539 genes")
expressed.id = which(apply(log(rpm+1,2),1,mean)!=0)
rpm = rpm[expressed.id,]

#----------------------------------------------------------------


#Table S1. ANOVA variation explaind by genotypes-----------------

options(contrasts=c("contr.sum", "contr.poly")) # set Type III ANOVA

v.list = c(); p.list = c()
for(i in 1:nrow(rpm)) { print(i)
  d = cbind(rpm[i,],primerID[,-c(1:12)])
  aov.res = try(aov(log(d[,1]+1,2)~Line+LeafLen+Bolting,data = d))
  
  if(class(aov.res)=="try-error") {v = NA} else {
    v = (drop1(aov.res,~.,test="F")$"Sum of Sq"[2])/(drop1(aov.res,~.,test="F")$RSS[2])
    p = drop1(aov.res,~.,test="F")$Pr[2]
  }
  v.list = c(v.list,v); p.list = c(p.list,p)
}

p.list = p.adjust(p.list,method="BH")

v.list = cbind(des[expressed.id,c(1:5,7)],v.list,p.list)
write.csv(v.list,"SatoAthRNA_varList190418.csv",row.names = F)
#----------------------------------------------------------------


#Table 2. GO enrichment analysis----------------------------------
library(GO.db) #load a R package

ng.mft = function(
  cgt, #output from ng.MakeContigGOidTable, [,1]:"locus", [,2]:"GOid"
  gn.test, #contig names for test
  alternative="greater"
){
  
  #cat(sprintf("%s\n", Sys.time()))
  
  gid.u = unique(cgt[,"GOid"])
  
  ft.in = matrix(0, nrow=length(gid.u), ncol=9)
  colnames(ft.in) = c("xtt", "xft", "xtf", "xff", "xnt", "xnf", "xtn", "xfn", "xnn")
  rownames(ft.in) = gid.u
  
  #               gn.test
  #             TRUE FALSE
  #Group  TRUE   xtt   xft   xnt
  #      FALSE   xtf   xff   xnf
  #              xtn   xfn   xnn
  #output = c("xtt", "xtn", "xnt", "xnn")
  
  ft.in[,"xnn"] = length(unique(cgt[, "locus"]))
  
  gn.pp.gid = table(cgt[, "GOid"])
  ft.in[names(gn.pp.gid), "xnt"] = gn.pp.gid
  ft.in[,"xnf"] = ft.in[,"xnn"] - ft.in[,"xnt"]
  
  ft.in[,"xtn"] = length(intersect(gn.test, unique(cgt[, "locus"])))
  ft.in[,"xfn"] = ft.in[,"xnn"] - ft.in[,"xtn"]
  
  gsea.test = cgt[is.element(cgt[,"locus"], gn.test), ]
  gn.test.gid = table(gsea.test[, "GOid"])
  ft.in[names(gn.test.gid), "xtt"] = gn.test.gid
  
  ft.in[,"xtf"] = ft.in[,"xtn"] - ft.in[,"xtt"]
  ft.in[,"xft"] = ft.in[,"xnt"] - ft.in[,"xtt"]
  ft.in[,"xff"] = ft.in[,"xnf"] - ft.in[,"xtf"]
  
  #cat(sprintf("%s\n", Sys.time()))
  
  #Fisher's exact test.  8? sec
  fr = rep(1, nrow(ft.in))
  dt = rep(1, nrow(ft.in))
  for(i in 1:nrow(ft.in)){
    start = Sys.time()
    if(ft.in[i,"xtn"] > 1 && ft.in[i,"xnt"] > 1){ 
      contable = matrix(ft.in[i, 1:4], ncol=2)
      tmp = fisher.test(contable, alternative = alternative)
      fr[i] = tmp$p.value
    } else {
    }
    end = Sys.time()
    dt[i] = end - start
  }
  
  out = cbind(fr, ft.in, dt)
  colnames(out) = c("p.value", colnames(ft.in), "time")
  rownames(out) = rownames(ft.in)
  
  #cat(sprintf("%s\n", Sys.time()))
  
  return(out)
  
}

# get GO terms -----------------------------
ng.GetGOTerms = function(GOid){
  
  out = NULL
  for(i in GOid){
    tmp = try(Term(i), silent=TRUE)
    if(class(tmp)=="try-error"){
      out = c(out, "NA")
    } else {
      out = c(out, Term(i))
    }
  }
  return(out)
}

# GO.slim--------------------
GO.slim = function(x){
  offspring = c(as.list(GOBPOFFSPRING),as.list(GOMFOFFSPRING),as.list(GOCCOFFSPRING))
  out = NULL
  for(i in 1:length(x)){
    if(sum(is.element(offspring[[x[i]]],x))==0){
      out = c(out, x[i])
    }
  }
  return(out)
}

# prep. GO test output table ------------------------------
ng.prepGOtestOutTable = function(r){
  
  adp = p.adjust(r[,"p.value"], method="BH")
  
  tmp.id = rownames(r)[adp<0.05]
  tmp.adp = adp[adp<0.05]
  tmp.description = ng.GetGOTerms(tmp.id)
  tmp.xnn = r[adp<0.05, c("xtt", "xtn", "xnt", "xnn")]
  
  out = cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  colnames(out) = c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  if(length(tmp.id)>0){
    return(out[is.element(tmp.id,GO.slim(tmp.id)),])
  }
}

GO.category = function(GO.list){
  BP = as.list(GOBPANCESTOR)
  MF = as.list(GOMFANCESTOR)
  CC = as.list(GOCCANCESTOR)
  GO.BP = GO.list[is.element(GO.list,names(BP))]
  GO.MF = GO.list[is.element(GO.list,names(MF))]
  GO.CC = GO.list[is.element(GO.list,names(CC))]
  return(list(BP=GO.BP, MF=GO.MF,CC=GO.CC))
}

#prepare a list of the top 5% genes showing a large among-accession variation
v.sum = read.csv("SatoAthRNA_varList190418.csv",header=T)
v.uplimit = quantile(na.omit(v.sum$v.list),0.95)
v.u95.list = v.sum[v.sum$v.list>v.uplimit,]
gl = v.u95.list$ID

#fisher tests
fisher.res = ng.mft(ulg, gl)
GOtable = ng.prepGOtestOutTable(fisher.res)
GOcat = GO.category(GOtable[,2])
write.csv(GOtable[GOcat$BP,],"SatoAthRNA_GOslimBP190418.csv",row.names = F)
#----------------------------------------------------------------


#Table S2. G x H test---------------------------------------------

options(contrasts=c("contr.sum", "contr.poly")) # set Type III ANOVA

#Le: L. erysimi; Ps, Pa: P. striolata and P. atra; Px: P. xylostella; Fo: F. occidentalis; Hole: Leaf holes made by Ps and Pa
p.list=c()
for(i in v.u95.list$ID) { print(i)
  d = cbind(rpm[i,],primerID)
  
  aov.main = aov(log(d[,1]+1,2)~Line+LeafLen+Bolting+Line+log(Le_sum+1),data = d)
  p1 = drop1(aov.main,~.,test="F")$Pr[5]
  aov.main = aov(log(d[,1]+1,2)~Line+LeafLen+Bolting+Line*log(Le_sum+1),data = d)
  p2 = drop1(aov.main,~.,test="F")$Pr[6]
  
  P_sum = d$Ps + d$Pa
  aov.int = aov(log(d[,1]+1,2)~Line+LeafLen+Bolting+Line+log(P_sum+1),data = d)
  p3 = drop1(aov.int,~.,test="F")$Pr[5]
  aov.int = aov(log(d[,1]+1,2)~Line+LeafLen+Bolting+Line*log(P_sum+1),data = d)
  p4 = drop1(aov.int,~.,test="F")$Pr[6]
  
  aov.int = aov(log(d[,1]+1,2)~Line+LeafLen+Bolting+Line+log(Px+1),data = d)
  p5 = drop1(aov.int,~.,test="F")$Pr[5]
  aov.int = aov(log(d[,1]+1,2)~Line+LeafLen+Bolting+Line*log(Px+1),data = d)
  p6 = drop1(aov.int,~.,test="F")$Pr[6]
  
  aov.int = aov(log(d[,1]+1,2)~Line+LeafLen+Bolting+Line+log(Fo+1),data = d)
  p7 = drop1(aov.int,~.,test="F")$Pr[5]
  aov.int = aov(log(d[,1]+1,2)~Line+LeafLen+Bolting+Line*log(Fo+1),data = d)
  p8 = drop1(aov.int,~.,test="F")$Pr[6]
  
  aov.int = aov(log(d[,1]+1,2)~Line+LeafLen+Bolting+Line+log(Hole+1),data = d)
  p9 = drop1(aov.int,~.,test="F")$Pr[5]
  aov.int = aov(log(d[,1]+1,2)~Line+LeafLen+Bolting+Line*log(Hole+1),data = d)
  p10 = drop1(aov.int,~.,test="F")$Pr[6]
  
  p.vec = c(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
  p.list = rbind(p.list, p.vec)
}

for(j in 1:ncol(p.list)) { #p.adjust by FDR
  p.list[,j] = p.adjust(p.list[,j],method="BH")
}

p.list = cbind(as.character(v.u95.list$ID),as.character(v.u95.list$Short.description),as.character(v.u95.list$gene.alias),v.u95.list$v.list,p.list)
colnames(p.list) = c("GeneID", "Description", "Alias", "v", "Le_main", "Le_int", "P_main", "P_int", "Px_main", "Px_int", "Fo_main", "Fo_int", "hole_main", "hole_int")
write.csv(p.list,"SatoAthRNA_GxH_plist190418.csv",row.names = F)
#----------------------------------------------------------------


#Figures-------------------------------------------------

#Figure 2. GSL variation
par(mfcol=c(4,2))

hist(v.sum$v.list,main=paste("upper 95%CI = ",quantile(na.omit(v.sum$v.list),0.95)), breaks=seq(from = 0, to = 0.85, by=0.025),xlab="Proportion of variation explained")
abline(v=quantile(na.omit(v.sum$v.list),0.95),lty=2)

gene.id = c("AT4G03050","AT4G03060","AT3G14210","AT5G26000","AT5G25980","AT1G54040","AT5G07690") #"AT4G03050" AOP3; "AT4G03060" AOP2; "AT3G14210" ESM1; "AT5G26000" TGG1; "AT5G25980" TGG2; "AT1G54040" ESP; "AT5G07690" MYB29

for(i in gene.id) {
  d = cbind(rpm[i,],primerID)
  lm.res = lm(log(d[,1]+1,2)~Line+Line*log(Le_sum+1,2),data = d)
  
  mean.rpm = cbind(aggregate(log(d[,1]+1,2)~Line,data=d,mean),aggregate(log(d[,1]+1,2)~Line,data=d,sd)[,2]/sqrt(aggregate(log(d[,1]+1,2)~Line,data=d,length)[,2]))
  mean.rpm = mean.rpm[order(mean.rpm[,2]),]
  
  mean.bar = barplot(mean.rpm[,2],ylim=c(0,max(mean.rpm[,2])+max(mean.rpm[,3])),names=mean.rpm[,1],las=2,ylab="log2(rpm + 1)",cex.names=1,main=i)
  arrows(mean.bar, mean.rpm[,2]+mean.rpm[,3], mean.bar, mean.rpm[,2]-mean.rpm[,3],code=3,length=0)
}

#Figure 3. AOP3-aphids, CYP81D11-beetles
par(mfcol=c(1,2))

#AOP3
gene.id = "AT4G03050" 
d = cbind(rpm[gene.id,],primerID)

name.l = "Ler-1"
plot(log(subset(d, Line==name.l)$Le_sum+1,2),jitter(log(subset(d, Line==name.l)[,1]+1,2)),cex=1.5,las=1,pch=16,ylim=c(0,12),xlab="log2(no. of aphids + 1)",ylab="log2(rpm+1)",main=name.l,col=grey(0.25,0.5))
name.l = "gl1-2"
points(log(subset(d, Line==name.l)$Le_sum+1,2),jitter(log(subset(d, Line==name.l)[,1]+1,2)),cex=1.5,las=1,pch=1,col=grey(0.5,0.5))

#CYP81D11
gene.id = "AT3G28740"
d = cbind(rpm[gene.id,],primerID)

plot(log(d$Ps+d$Pa + 1, 2), jitter(lm(log(d[,1]+1,2)~d$Line)$resid), cex=1.2, las=1, pch=16, col=grey(0.5,0.3), main=gene.id, ylab="resid. log2(rpm + 1)", xlab="log2(No. of Phyllotreta beetles + 1)")
reg = lm(lm(log(d[,1]+1,2)~d$Line)$resid~log(d$Ps+d$Pa + 1, 2))
abline(a=reg$coefficients[1], b=reg$coefficients[2])

#-----------------------------------------------------------------------

#Figure S1: Comparison with Chan et al. (2010) GSL data-----------------

d.gls = read.csv("Chan2010GWAS_averageGLS.csv", header=T)

#AOP3, AT4G03050; 3 OH-propyl
gene.id = "AT4G03050"
d = cbind(rpm[gene.id,],primerID)

mean.rpm = aggregate(log(d[,1]+1,2)~Line,data=d,mean)
mean.rpm = mean.rpm[which(mean.rpm[,1]!="gl1-1"&mean.rpm[,1]!="gl1-2"),] #excl. gl1 mutants

gls.line=c()
for(i in mean.rpm[,1]) {
  gls.line = rbind(gls.line, subset(d.gls, Accession==i))
}

cor.test(log(gls.line[,2]+10^-2,2),mean.rpm[,2])
plot(log(gls.line[,2]+10^-2,2),mean.rpm[,2], las=1, pch=16)


#AOP2, AT4G03060: Alkenyl
gene.id = "AT4G03060"
d = cbind(rpm[gene.id,],primerID)

mean.rpm = aggregate(log(d[,1]+1,2)~Line,data=d,mean)
mean.rpm = mean.rpm[which(mean.rpm[,1]!="gl1-1"&mean.rpm[,1]!="gl1-2"),] #excl. gl1 mutants

cor.test(log(gls.line[,8]+10^-2,2),mean.rpm[,2])
plot(log(gls.line[,8]+10^-2,2),mean.rpm[,2],las=1)

#----------------------------------------------------------------

