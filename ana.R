#=====================================================
#ana.R
#=====================================================
library(ggpubr)
library(ggrepel) 
library(RColorBrewer)
library(forestplot)
library(grid)
library(dplyr)

setwd("../files/")
dis1<-read.csv("../files/supplemental_model1.csv") %>%
  mutate(Estimate=beta_dis, P=padj_dis,thres=padj_dis<0.05) %>%
  dplyr::select(SomaSeqid, Target, Estimate, P, thres) %>%
  arrange(P)
dis1$label<-""
dis1$label[1:10]<-as.character(dis1$Target[1:10])
dis2<-read.csv("../files/supplemental_model2.csv") %>%
  mutate(Estimate=beta_dis, P=padj_dis,thres=padj_dis<0.05) %>%
  dplyr::select(SomaSeqid, Target, Estimate, P, thres) %>%
  arrange(P)
dis2$label<-""
dis2$label[1:10]<-as.character(dis2$Target[1:10])
dis3<-read.csv("../files/supplemental_model3.csv") %>%
  mutate(Estimate=beta_dis, P=padj_dis,thres=padj_dis<0.05) %>%
  dplyr::select(SomaSeqid, Target, Estimate, P, thres) %>%
  arrange(P)
dis3$label<-""
dis3$label[1:10]<-as.character(dis3$Target[1:10])

rep1<-read.csv("../files/supplemental_model1.csv") %>%
  mutate(Estimate=beta_rep, P=padj_rep,thres=padj_rep<0.05) %>%
  dplyr::select(SomaSeqid, Target, Estimate, P, thres)%>%
  arrange(P)
rep1$label<-""
rep1$label[1:10]<-as.character(rep1$Target[1:10])
rep2<-read.csv("../files/supplemental_model2.csv") %>%
  mutate(Estimate=beta_rep, P=padj_rep,thres=padj_rep<0.05) %>%
  dplyr::select(SomaSeqid, Target, Estimate, P, thres)%>%
  arrange(P)
rep2$label<-""
rep2$label[1:10]<-as.character(rep2$Target[1:10])
rep3<-read.csv(../files/supplemental_model3.csv") %>%
  mutate(Estimate=beta_rep, P=padj_rep,thres=padj_rep<0.05) %>%
  dplyr::select(SomaSeqid, Target, Estimate, P, thres)%>%
  arrange(P)
rep3$label<-""
rep3$label[1:10]<-as.character(rep3$Target[1:10])

#========2022-05-11 should use p_meta instead of padj_meta
meta1<-read.csv("../files/supplemental_model1.csv") %>%
  #mutate(Estimate=beta_meta, P=padj_meta,thres=padj_meta<(0.05/4301)) %>%
  mutate(Estimate=beta_meta, P=p_meta,thres=p_meta<(0.05/4301)) %>%
  dplyr::select(SomaSeqid, Target, Estimate, P, thres)%>%
  arrange(P)
meta1$label<-""
meta1$label[1:10]<-as.character(meta1$Target[1:10])
meta2<-read.csv("../files/supplemental_model2.csv") %>%
  #mutate(Estimate=beta_meta, P=padj_meta,thres=padj_meta<(0.05/4301)) %>%
  mutate(Estimate=beta_meta, P=p_meta,thres=p_meta<(0.05/4301)) %>%
  dplyr::select(SomaSeqid, Target, Estimate, P, thres)%>%
  arrange(P)
meta2$label<-""
meta2$label[1:10]<-as.character(meta2$Target[1:10])
meta3<-read.csv("../files//supplemental_model3.csv") %>%
  #mutate(Estimate=beta_meta, P=padj_meta,thres=padj_meta<(0.05/4301)) %>%
  mutate(Estimate=beta_meta, P=p_meta,thres=p_meta<(0.05/4301)) %>%
  dplyr::select(SomaSeqid, Target, Estimate, P, thres)%>%
  arrange(P)
meta3$label<-""
meta3$label[1:10]<-as.character(meta3$Target[1:10])

savem1.out<-list()
for (model in 1:3) {
  savem1.out[[model]]<-vector("list",dim(dis1)[1])
  savem1.out[[model]]<-get(paste0("dis",model))
  #assign(savem1.out[[model]],paste0("dis",model))
}
names(savem1.out)<-c("outcome1","outcome2","outcome3") 

savem2.out<-list()
for (model in 1:3) {
  savem2.out[[model]]<-vector("list",dim(rep1)[1])
  savem2.out[[model]]<-get(paste0("rep",model))
}
names(savem2.out)<-c("outcome1","outcome2","outcome3") 

savem3.out<-list()
for (model in 1:3) {
  savem3.out[[model]]<-vector("list",dim(meta1)[1])
  savem3.out[[model]]<-get(paste0("meta",model))
}
names(savem3.out)<-c("outcome1","outcome2","outcome3") 


save.all<-list()
save.all[[1]]<-savem1.out
save.all[[2]]<-savem2.out
save.all[[3]]<-savem3.out

names(save.all)<-c("Discovery","Replication","Meta")

analysis<-c("Discovery","Replication","Meta")
vol = vector("list",3)
names(vol) = analysis
for (ana in 1:3) {
  vol[[ana]]<-vector("list",3)
  names(vol[[ana]])<-c("outcome1","outcome2","outcome3")
  for (i in 1:3) {
    out=save.all[[ana]][[i]]
    ##
    vol[[ana]][[i]] = ggplot(out) + 
      geom_point(aes(x=Estimate, y=-log10(P), colour=thres)) +
      scale_color_manual(values = c('gray','blue')) + 
      geom_text_repel(aes(x = Estimate, y =-log10(P), label = label),size = 2.5, box.padding = unit(0.5, "lines"), max.overlaps = Inf) +
      theme_classic() + xlab('Effect size') + ylab('-log10(P-value)') + theme(legend.position = "none", text = element_text(size=10)) 
  }
}


ggarrange(vol$Discovery$outcome1 + ggtitle(''), 
          vol$Discovery$outcome2 + ggtitle(''), 
          vol$Discovery$outcome3 + ggtitle(''), 
          nrow=3, ncol=1)
ggsave(paste0("volcano_discovery",Sys.Date(),".png"),device = "png", width = 4, height = 11)

ggarrange(vol$Replication$outcome1 + ggtitle(''), 
          vol$Replication$outcome2 + ggtitle(''), 
          vol$Replication$outcome3 + ggtitle(''), 
          nrow=3, ncol=1)
ggsave(paste0("volcano_replication",Sys.Date(),".png"),device = "png", width = 4, height = 11)

ggarrange(vol$Meta$outcome1 + ggtitle(''), 
          vol$Meta$outcome2 + ggtitle(''), 
          vol$Meta$outcome3 + ggtitle(''), 
          nrow=3, ncol=1)
ggsave(paste0("volcano_meta",Sys.Date(),".png"),device = "png", width = 4, height = 11)
