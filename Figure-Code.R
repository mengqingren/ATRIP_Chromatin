setwd("/mnt/sda/Data_ZhuQian/PLA/PLA")
library(ggplot2)
library(ggbeeswarm)
colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c','#7bc7cd','#5d84a4','#4B4B5D',"#EC7232")

colorlist2 = c("#ab586c","#efe9b6","#9ebeb9","#f5b366","#88758f","#c8b0a6","#c29367","#fbbb8d","#b48b87","#85a283",
               "#e08988","#c7828b","#c1c1ca","#d3b171","#d7c4c9","#65a9a7","#6f7d9b","#9b8651","#6f8d93","#6e7d9c","#706b70",
               "#716c72","#d47e4d","#995181","#ad9eba","#e08593","#876f69","#c89e6b","#ecdf98","#82b1c3","#bf9e97","#7ab884",
               "#8b5881","#97b479","#cfac68","#eabfa8","#e4e09f","#64a0ab","#809163","#ad6c87","#629260","#cb905e","#b2613f",
               "#b0d8dc","#548550","#d38448","#b65b66","#618d58","#7ba487","#f3cea8","#7ba38b","#ba5b63")
colorlist3 = c("#ab586c","#9ebeb9","#f5b366","#88758f","#c8b0a6","#c29367")
############ EdU-CDC45 ###########
library(readxl)
dat <- read_xlsx("EdU-CDC45/CTR shATRIP shHP1B shSUV420H2.xlsx")
dat2 <- dat %>% gather(group,Count) %>%
  mutate(group=factor(group,levels=colnames(dat)))%>%
  filter(Count<60)
dat3 <- dat2 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "EdU-CDC45/CTR-shATRIP-shHP1B-shSUV420H2.Wilcox.csv",row.names = F)
p1 <- ggplot(dat2, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  #ylim(0,100)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=34, yend=34,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=36,label="n.s.",size=4)+
  annotate("segment", x=2, xend=3, y=37, yend=37,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=37.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=40, yend=40,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=40.5,label="****",size=4)+
  annotate("segment", x=2, xend=5, y=43, yend=43,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=43.5,label="***",size=4)
ggsave(p1,filename="EdU-CDC45/CTR-shATRIP-shHP1B-shSUV420H2.pdf",height = 2,width = 4)


library(readxl)
dat <- read_xlsx("EdU-CDC45/EdU-CDC45 (shATRIP rescue).xlsx")
dat2 <- dat %>% gather(group,Count) %>%
  mutate(group=factor(group,levels=colnames(dat)))
dat3 <- dat2 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "EdU-CDC45/EdU-CDC45-shATRIP-rescue.Wilcox.csv",row.names = F)
p2 <- ggplot(dat2, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,80)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=65, yend=65,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=67,label="n.s.",size=3)+
  annotate("segment", x=2, xend=3, y=68, yend=68,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=68.5,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=68, yend=68,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=68.5,label="****",size=4)+
  annotate("segment", x=4, xend=5, y=68, yend=68,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=68.5,label="****",size=4)+
  annotate("segment", x=3, xend=5, y=73, yend=73,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=75,label="n.s.",size=3)+
  annotate("segment", x=2, xend=5, y=78, yend=78,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=78.5,label="****",size=4)
ggsave(p2,filename="EdU-CDC45/EdU-CDC45-shATRIP-rescue.pdf",height = 2,width = 4)

############ EdU-MCM2 ###########
dat <- read_xlsx("EdU-MCM2/CTR shATRIP shHP1B shSUV420H2.xlsx")
dat2 <- dat %>% gather(group,Count) %>%
  mutate(group=factor(group,levels=colnames(dat))) %>%
  filter(Count<40)
dat3 <- dat2 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "EdU-MCM2/CTR-shATRIP-shHP1B-shSUV420H2.Wilcox.csv",row.names = F)
p1 <- ggplot(dat2, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  #ylim(0,100)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=34, yend=34,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=36,label="n.s.",size=4)+
  annotate("segment", x=2, xend=3, y=37, yend=37,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=37.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=40, yend=40,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=40.5,label="****",size=4)+
  annotate("segment", x=2, xend=5, y=43, yend=43,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=43.5,label="****",size=4)
ggsave(p1,filename="EdU-MCM2/CTR-shATRIP-shHP1B-shSUV420H2.pdf",height = 2,width = 4)

library(readxl)
dat <- read_xlsx("EdU-MCM2/EdU-MCM2 (shATRIP rescue).xlsx")
dat2 <- dat %>% gather(group,Count) %>%
  mutate(group=factor(group,levels=colnames(dat)))
dat3 <- dat2 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "EdU-MCM2/EdU-MCM2-shATRIP-rescue.Wilcox.csv",row.names = F)
p2 <- ggplot(dat2, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,80)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=65, yend=65,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=67,label="n.s.",size=3)+
  annotate("segment", x=2, xend=3, y=68, yend=68,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=68.5,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=68, yend=68,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=68.5,label="****",size=4)+
  annotate("segment", x=4, xend=5, y=68, yend=68,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=68.5,label="****",size=4)+
  annotate("segment", x=3, xend=5, y=73, yend=73,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=75,label="n.s.",size=3)+
  annotate("segment", x=2, xend=5, y=78, yend=78,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=78.5,label="****",size=4)
ggsave(p2,filename="EdU-MCM2/EdU-MCM2-shATRIP-rescue.pdf",height = 2,width = 4)

############ F2.B ###########
dat1 <- read_xlsx("F2.B/EdU&ATRIP.xlsx") %>% gather(group,Count) %>% mutate(G="PLA:EdU-ATRIP")
dat2 <- read_xlsx("F2.B/EdU&HP1B.xlsx") %>% gather(group,Count) %>% mutate(G="PLA:EdU-HP1B")
dat3 <- read_xlsx("F2.B/EdU&SUV420H2.xlsx") %>% gather(group,Count) %>% mutate(G="PLA:EdU-SUV420H2")
dat4 <- read_xlsx("F2.B/EdU&H4K20me3.xlsx") %>% gather(group,Count) %>% mutate(G="PLA:EdU-H4K20me3")
dat <- rbind(dat1,dat2) %>% rbind(dat3) %>% rbind(dat4)
write.csv(dat,"F2.B/EdU-ATRIP-HP1B-SUV420H2-H4K20me3.csv",row.names = F,quote = F)

p4.1 <- ggplot(dat %>% filter(G == "PLA:EdU-ATRIP"), 
       aes(x = group, y = Count, color = group)) +
  geom_violin(aes(fill=group),alpha=0.6)+
  #geom_boxplot(width=0.01,fill="grey")+
  #geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.2,jitter.height = 0.000001),size=1.5,alpha=0.4)+
  geom_quasirandom(size = 1.5, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(aes(color=group),shape=0,fun = "mean",size=2,color="black",geom = "point")+
  stat_summary(fun.data = "mean_se",geom = "errorbar",color="black",width = 0)+
  scale_color_manual("",values = colorlist1) +
  scale_fill_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "PLA foci per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,25)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  geom_signif(comparisons = list(c("CTR","HU")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("PLA\nEdU-ATRIP")
p4.2 <- ggplot(dat %>% filter(G == "PLA:EdU-HP1B"), 
               aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "PLA foci per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  ylim(0,42)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  geom_signif(comparisons = list(c("CTR","HU")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("PLA\nEdU-HP1B")
p4.3 <- ggplot(dat %>% filter(G == "PLA:EdU-SUV420H2"), 
               aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "PLA foci per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  ylim(0,38)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  geom_signif(comparisons = list(c("CTR","HU")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("PLA\nEdU-SUV4-20H2")
p4.4 <- ggplot(dat %>% filter(G == "PLA:EdU-H4K20me3"), 
               aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "PLA foci per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  ylim(0,38)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  geom_signif(comparisons = list(c("CTR","HU")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("PLA\nEdU-H4K20me3")
library(cowplot)
p4 <- plot_grid(p4.1,p4.2,p4.3,p4.4,nrow = 1)
ggsave(p4,filename="F2.B/EdU-ATRIP-HP1B-SUV420H2-H4K20me3.pdf",height = 2.4,width = 10)

############ F3.C ###########
dat1 <- read_excel("F3.C/ATRIP&HP1B.xlsx") %>% dplyr::select(-...1) %>%
  gather(group,Count) %>% na.omit()
p5.2 <- ggplot(dat1,aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "PLA foci per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,30)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  geom_signif(comparisons = list(c("CTR","HU")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("PLA\nATRIP-HP1B")

dat2 <- read_excel("F3.C/EdU&ATRIP.xlsx") %>% dplyr::select(-...1) %>%
  gather(group,Count) %>% na.omit()
p5.1 <- ggplot(dat2,aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "PLA foci per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,30)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  geom_signif(comparisons = list(c("CTR","HU")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("PLA\nEDU-ATRIP")

dat3 <- read_excel("F3.C/ATRIP&SUV420H2.xlsx") %>% dplyr::select(-...1) %>%
  gather(group,Count) %>% na.omit()
p5.3 <- ggplot(dat2,aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "PLA foci per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,30)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  geom_signif(comparisons = list(c("CTR","HU")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("PLA\nATRIP-SUV4-20H2")

library(cowplot)
p5 <- plot_grid(p5.1,p5.2,p5.3,nrow = 1)
ggsave(p5,filename="F3.C/PLA-EdU-ATRIP-HP1B-SUV420H2.pdf",height = 2.4,width = 7.5)

############ F3.D ###########
dat1 <- read_xlsx("F3.D/EdU&ATRIP.xlsx") %>%
  magrittr::set_colnames(c("NC HU","R-shATRIP-Vector","R-shATRIP-WT","R-shATRIP-Del")) %>%
  gather(group,Count) %>% na.omit() %>%
  mutate(group=factor(group,levels=c("NC HU","R-shATRIP-Vector","R-shATRIP-WT","R-shATRIP-Del")))
dat3 <- dat1 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "F3.D/EdU&ATRIP.Wilcox.csv",row.names = F)
p6 <- ggplot(dat1,aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist3) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5,colour = "black"),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,45)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  #geom_signif(comparisons = list(c("shATRIP","R-shATRIP-Vector"),c("shATRIP","R-shATRIP-WT"),c("shATRIP","R-shATRIP-Del"),c("R-shATRIP-WT","R-shATRIP-Vector"),c("R-shATRIP-WT","R-shATRIP-Del")),map_signif_level = T,test =wilcox.test)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("PLA\nEdU-ATRIP")+
  annotate("segment", x=1, xend=2, y=16, yend=16,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=17,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=40, yend=40,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=41,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=40, yend=40,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=42,label="n.s.",size=3)+
  annotate("segment", x=2, xend=4, y=43, yend=43,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=43.5,label="***",size=4)
ggsave(p6,filename="F3.D/EdU-ATRIP.pdf",height = 2.4,width = 4)

############ F3.J ###########
dat1 <- read_xlsx("F3.J/EdU-H4K20me3 shATRIP rescue.xlsx") %>% gather(group,Count) %>% na.omit() %>%
  mutate(group=factor(group,levels=c("NC CTR","NC HU","shATRIP-Vector HU","shATRIP-WT HU","shATRIP-Del HU")))
dat3 <- dat1 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "F3.J/EdU-H4K20me3-shATRIPrescue.Wilcox.csv",row.names = F)

p1 <- ggplot(dat1, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,80)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=40, yend=40,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=41,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=40, yend=40,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=41,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=64, yend=64,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=64.5,label="****",size=4)+
  annotate("segment", x=4, xend=5, y=64, yend=64,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=64.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=72, yend=72,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=72.5,label="****",size=4)+
  annotate("segment", x=3, xend=5, y=68, yend=68,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=68.5,label="****",size=4)+
  annotate("segment", x=2, xend=5, y=76, yend=76,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=76.5,label="****",size=4)+
  ggtitle("PLA\nEdU-H4K20me3")
ggsave(p1,filename="F3.J/EdU-H4K20me3-shATRIPrescue.pdf",height = 2.4,width = 4.1)

dat1 <- read_xlsx("F3.J/EdU-HP1B shATRIP rescue.xlsx")%>% gather(group,Count) %>% na.omit() %>%
  mutate(group=factor(group,levels=c("NC CTR","NC HU","shATRIP-Vector HU","shATRIP-WT HU","shATRIP-Del HU")))
dat3 <- dat1 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "F3.J/EdU-HP1B-shATRIPrescue.Wilcox.csv",row.names = F)

p1 <- ggplot(dat1, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,80)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=64, yend=64,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=64.5,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=64, yend=64,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=64.5,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=64, yend=64,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=64.5,label="****",size=4)+
  annotate("segment", x=4, xend=5, y=64, yend=64,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=64.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=72, yend=72,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=74,label="n.s.",size=3)+
  annotate("segment", x=3, xend=5, y=68, yend=68,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=70,label="n.s.",size=3)+
  annotate("segment", x=2, xend=5, y=76, yend=76,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=76.5,label="****",size=4)+
  ggtitle("PLA\nEdU-HP1B")
ggsave(p1,filename="F3.J/EdU-HP1B-shATRIPrescue.pdf",height = 2.4,width = 4.1)

dat1 <- read_xlsx("F3.J/EdU-SUV420H2 shATRIP rescue.xlsx") %>% gather(group,Count) %>% na.omit() %>%
  mutate(group=factor(group,levels=c("NC CTR","NC HU","shATRIP-Vector HU","shATRIP-WT HU","shATRIP-Del HU")))
dat3 <- dat1 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "F3.J/EdU-SUV420H2-shATRIPrescue.Wilcox.csv",row.names = F)

p1 <- ggplot(dat1, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,85)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=43, yend=43,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=43.5,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=43, yend=43,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=43.5,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=72, yend=72,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=72.5,label="****",size=4)+
  annotate("segment", x=4, xend=5, y=72, yend=72,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=72.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=80, yend=80,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=80.5,label="*",size=4)+
  annotate("segment", x=3, xend=5, y=76, yend=76,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=76.5,label="*",size=4)+
  annotate("segment", x=2, xend=5, y=84, yend=84,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=84.5,label="****",size=4)+
  ggtitle("PLA\nEdU-SUV4-20H2")
ggsave(p1,filename="F3.J/EdU-SUV420H2-shATRIPrescue.pdf",height = 2.4,width = 4.1)

############ F3.K ###########
dat1 <- read_xlsx("F3.K/EdU-H4K20me3(shHP1B rescue).xlsx") %>% gather(group,Count) %>% na.omit() %>%
  mutate(group=factor(group,levels=c("NC CTR","NC HU","shHP1β-Vector","shHP1β-WT","shHP1β-ΔCSD")))
dat3 <- dat1 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "F3.K/EdU-H4K20me3-shHP1Brescue.Wilcox.csv",row.names = F)

p1 <- ggplot(dat1, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,60)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=50, yend=50,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=51,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=50, yend=50,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=51,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=40, yend=40,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=40.5,label="****",size=4)+
  annotate("segment", x=4, xend=5, y=40, yend=40,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=40.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=44, yend=44,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=46,label="n.s.",size=3)+
  annotate("segment", x=3, xend=5, y=48, yend=48,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=50,label="n.s.",size=3)+
  annotate("segment", x=2, xend=5, y=54, yend=54,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=54.5,label="****",size=4)+
  ggtitle("PLA\nEdU-H4K20me3")
ggsave(p1,filename="F3.K/EdU-H4K20me3-shHP1Brescue.pdf",height = 2.4,width = 4.1)

dat1 <- read_xlsx("F3.K/EdU-HP1B (shHP1B rescue).xlsx") %>% gather(group,Count) %>% na.omit() %>%
  mutate(group=factor(group,levels=c("NC CTR","NC HU","shHP1β-Vector","shHP1β-WT","shHP1β-ΔCSD")))
dat3 <- dat1 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "F3.K/EdU-HP1B-shHP1Brescue.Wilcox.csv",row.names = F)

p1 <- ggplot(dat1, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,55)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=41, yend=41,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=41.5,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=41, yend=41,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=41.5,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=36, yend=36,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=38,label="n.s.",size=3)+
  annotate("segment", x=4, xend=5, y=36, yend=36,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=36.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=44, yend=44,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=44.5,label="**",size=4)+
  annotate("segment", x=3, xend=5, y=48, yend=48,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=50,label="n.s.",size=3)+
  annotate("segment", x=2, xend=5, y=52, yend=52,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=52.5,label="****",size=4)+
  ggtitle("PLA\nEdU-HP1B")
ggsave(p1,filename="F3.K/EdU-HP1B-shHP1Brescue.pdf",height = 2.4,width = 4.1)

dat1 <- read_xlsx("F3.K/EdU-SUV420H2(shHP1B rescue).xlsx") %>% gather(group,Count) %>% na.omit() %>%
  mutate(group=factor(group,levels=c("NC CTR","NC HU","shHP1β-Vector","shHP1β-WT","shHP1β-ΔCSD")))
dat3 <- dat1 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "F3.K/EdU-SUV420H2-shHP1Brescue.Wilcox.csv",row.names = F)

p1 <- ggplot(dat1, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,52)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=37, yend=37,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=38,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=37, yend=37,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=38,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=40, yend=40,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=40.5,label="****",size=4)+
  annotate("segment", x=4, xend=5, y=40, yend=40,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=40.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=43, yend=43,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=45,label="n.s.",size=3)+
  annotate("segment", x=3, xend=5, y=46, yend=46,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=46.5,label="***",size=4)+
  annotate("segment", x=2, xend=5, y=49, yend=49,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=49.5,label="****",size=4)+
  ggtitle("PLA\nEdU-SUV4-20H2")
ggsave(p1,filename="F3.K/EdU-SUV420H2-shHP1Brescue.pdf",height = 2.4,width = 4.1)

############ F3.L ###########
dat1 <- read_xlsx("F3.L/EdU-H4K20me3(shSUV420H2rescue).xlsx") %>% 
  gather(group,Count) %>% na.omit() %>%
  mutate(group=factor(group,levels=c("NC CTR","NC HU","shSUV420-Vector","shSUV420H2-WT","shSUV420H2-Δclamp","shSUV420H2-ΔSET")))
dat3 <- dat1 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "F3.L/EdU-H4K20me3-shSUV420H2rescue.Wilcox.csv",row.names = F)

p1 <- ggplot(dat1, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,75)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=58, yend=58,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=58.5,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=58, yend=58,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=58.5,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=50, yend=50,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=50.5,label="****",size=4)+
  annotate("segment", x=4, xend=5, y=50, yend=50,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=50.5,label="****",size=4)+
  annotate("segment", x=5, xend=6, y=20, yend=20,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=5.5,y=20.5,label="*",size=4)+
  annotate("segment", x=2, xend=4, y=62, yend=62,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=62.5,label="***",size=4)+
  annotate("segment", x=2, xend=5, y=66, yend=66,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=66.5,label="****",size=4)+
  annotate("segment", x=2, xend=6, y=70, yend=70,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=70.5,label="****",size=4)+
  annotate("segment", x=3, xend=5, y=54, yend=54,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=56,label="n.s.",size=3)+
  annotate("segment", x=3, xend=6, y=58, yend=58,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=58.5,label="***",size=4)+
  annotate("segment", x=4, xend=6, y=62, yend=62,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=5,y=62.5,label="****",size=4)+
  ggtitle("PLA\nEdU-H4K20me3")
ggsave(p1,filename="F3.L/EdU-H4K20me3-shSUV420H2rescue.pdf",height = 2.4,width = 4.4)


dat1 <- read_xlsx("F3.L/EdU-HP1β（shSUV420H2 rescue).xlsx") %>% 
  gather(group,Count) %>% na.omit() %>%
  mutate(group=factor(group,levels=c("NC CTR","NC HU","shSUV420-Vector","shSUV420H2-WT","shSUV420H2-Δclamp","shSUV420H2-ΔSET")))
dat3 <- dat1 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "F3.L/EdU-HP1β-shSUV420H2rescue.Wilcox.csv",row.names = F)

p1 <- ggplot(dat1, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,60)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=40, yend=40,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=40.5,label="****",size=4)+
  ggtitle("PLA\nEdU-HP1β")
ggsave(p1,filename="F3.L/EdU-HP1β-shSUV420H2rescue.pdf",height = 2.4,width = 4.4)


dat1 <- read_xlsx("F3.L/EdU-SUV420H2(shSUV420H2 rescue).xlsx") %>% 
  gather(group,Count) %>% na.omit() %>%
  mutate(group=factor(group,levels=c("NC CTR","NC HU","shSUV420-Vector","shSUV420H2-WT","shSUV420H2-Δclamp","shSUV420H2-ΔSET")))
dat3 <- dat1 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "F3.L/EdU-SUV420H2-shSUV420H2rescue.Wilcox.csv",row.names = F)

p1 <- ggplot(dat1, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,80)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=65, yend=65,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=65.5,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=65, yend=65,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=65.5,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=50, yend=50,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=50.5,label="****",size=4)+
  annotate("segment", x=4, xend=5, y=50, yend=50,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=50.5,label="****",size=4)+
  annotate("segment", x=5, xend=6, y=48, yend=48,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=5.5,y=48.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=69, yend=69,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=69.5,label="****",size=4)+
  annotate("segment", x=2, xend=5, y=73, yend=73,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=73.5,label="****",size=4)+
  annotate("segment", x=2, xend=6, y=77, yend=77,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=77.5,label="****",size=4)+
  annotate("segment", x=3, xend=5, y=54, yend=54,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=56,label="n.s.",size=3)+
  annotate("segment", x=3, xend=6, y=58, yend=58,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=58.5,label="***",size=4)+
  annotate("segment", x=4, xend=6, y=62, yend=62,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=5,y=64,label="n.s.",size=3)+
  ggtitle("PLA\nEdU-SUV4-20H2")
ggsave(p1,filename="F3.L/EdU-SUV420H2-shSUV420H2rescue.pdf",height = 2.4,width = 4.4)

################ -------- #####################
setwd("/mnt/sda/Data_ZhuQian/Figure-LinePlot")
############ ATRIP-rescue => EdU-ATRIP ###########
colors19 = c("#b42e20","#a349a4","#313c63")
names(colors19) = c("H3","EdU","ATRIP")
dat <- read_xlsx("ATRIP-rescue/EdU-ATRIP-H3亮度统计（shATRIP rescue）/NC CTR/EdU-ATRIP NC CTR.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:ATRIP)
library(ggalt)

p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","ATRIP"))),
       aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="ATRIP-rescue/EdU-ATRIP-H3亮度统计（shATRIP rescue）/NC-CTR.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-ATRIP-H3亮度统计（shATRIP rescue）/NC HU/NC HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:ATRIP)
library(ggalt)

p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","ATRIP"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p2,filename="ATRIP-rescue/EdU-ATRIP-H3亮度统计（shATRIP rescue）/NC-HU.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-ATRIP-H3亮度统计（shATRIP rescue）/shATRIP-vector HU/shATRIP-vector HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:ATRIP)
library(ggalt)

p3 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","ATRIP"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p3,filename="ATRIP-rescue/EdU-ATRIP-H3亮度统计（shATRIP rescue）/shATRIP-vector-HU.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-ATRIP-H3亮度统计（shATRIP rescue）/shATRIP-Del HU/shATRIP-Del HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:ATRIP)
library(ggalt)

p3 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","ATRIP"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p3,filename="ATRIP-rescue/EdU-ATRIP-H3亮度统计（shATRIP rescue）/shATRIP-Del-HU.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-ATRIP-H3亮度统计（shATRIP rescue）/shATRIP-WT HU/shATRIP-WT HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:ATRIP)
library(ggalt)
p5 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","ATRIP"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p5,filename="ATRIP-rescue/EdU-ATRIP-H3亮度统计（shATRIP rescue）/shATRIP-WT-HU.pdf",height = 2,width = 4)

############ ATRIP-rescue => EdU-H4K20me3 ###########
colors19 = c("#b42e20","#a349a4","#5d84a4")
names(colors19) = c("H3","EdU","H4K20me3")
dat <- read_xlsx("ATRIP-rescue/EdU-H4K20me3-H3亮度统计（shATRIP rescue）/NC CTR/NC CTR.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:H4K20me3)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","H4K20me3"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="ATRIP-rescue/EdU-H4K20me3-H3亮度统计（shATRIP rescue）/NC-CTR.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-H4K20me3-H3亮度统计（shATRIP rescue）/NC HU/NC HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:H4K20me3)
library(ggalt)
p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","H4K20me3"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p2,filename="ATRIP-rescue/EdU-H4K20me3-H3亮度统计（shATRIP rescue）/NC-HU.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-H4K20me3-H3亮度统计（shATRIP rescue）/shATRIP-Del HU/shATRIP-Del HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:H4K20me3)
library(ggalt)
p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","H4K20me3"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p2,filename="ATRIP-rescue/EdU-H4K20me3-H3亮度统计（shATRIP rescue）/shATRIP-Del-HU.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-H4K20me3-H3亮度统计（shATRIP rescue）/shATRIP-vector HU/shATRIP-vector HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:H4K20me3)
library(ggalt)
p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","H4K20me3"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p2,filename="ATRIP-rescue/EdU-H4K20me3-H3亮度统计（shATRIP rescue）/shATRIP-vector-HU.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-H4K20me3-H3亮度统计（shATRIP rescue）/shATRIP-WT HU/shATRIP-WT HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:H4K20me3)
library(ggalt)
p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","H4K20me3"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p2,filename="ATRIP-rescue/EdU-H4K20me3-H3亮度统计（shATRIP rescue）/shATRIP-WT-HU.pdf",height = 2,width = 4)

############ ATRIP-rescue => EdU-HP1B ###########
colors19 = c("#b42e20","#a349a4","#7bc7cd")
names(colors19) = c("H3","EdU","HP1β")
dat <- read_xlsx("ATRIP-rescue/EdU-HP1B-H3亮度统计（shATRIP rescue）/NC CTR/NC CTR.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:HP1β)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","HP1β"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="ATRIP-rescue/EdU-HP1B-H3亮度统计（shATRIP rescue）/NC-CTR.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-HP1B-H3亮度统计（shATRIP rescue）/NC HU/NC HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:HP1β)
library(ggalt)
p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","HP1β"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p2,filename="ATRIP-rescue/EdU-HP1B-H3亮度统计（shATRIP rescue）/NC-HU.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-HP1B-H3亮度统计（shATRIP rescue）/shATRIP-Del HU/shATRIP-DEL HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:HP1β)
library(ggalt)
p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","HP1β"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p2,filename="ATRIP-rescue/EdU-HP1B-H3亮度统计（shATRIP rescue）/shATRIP-Del-HU.pdf",height = 2,width = 4)


dat <- read_xlsx("ATRIP-rescue/EdU-HP1B-H3亮度统计（shATRIP rescue）/shATRIP-vector HU/shATRIP-vector HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:HP1β)
library(ggalt)
p3 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","HP1β"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p3,filename="ATRIP-rescue/EdU-HP1B-H3亮度统计（shATRIP rescue）/shATRIP-vector-HU.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-HP1B-H3亮度统计（shATRIP rescue）/shATRIP-WT HU/shATRIP-WT HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:HP1β)
library(ggalt)
p3 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","HP1β"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p3,filename="ATRIP-rescue/EdU-HP1B-H3亮度统计（shATRIP rescue）/shATRIP-WT-HU.pdf",height = 2,width = 4)


############ ATRIP-rescue => EdU-SUV420H2 ###########
colors19 = c("#b42e20","#a349a4","#ebc03e")
names(colors19) = c("H3","EdU","SUV4-20H2")
dat <- read_xlsx("ATRIP-rescue/EdU-SUV420H2-H3亮度统计（shATRIP rescue）/NC CTR/NC CTR.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:`SUV4-20H2`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","SUV4-20H2"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="ATRIP-rescue/EdU-SUV420H2-H3亮度统计（shATRIP rescue）/NC-CTR.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-SUV420H2-H3亮度统计（shATRIP rescue）/NC HU/NC HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:`SUV4-20H2`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","SUV4-20H2"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="ATRIP-rescue/EdU-SUV420H2-H3亮度统计（shATRIP rescue）/NC-HU.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-SUV420H2-H3亮度统计（shATRIP rescue）/shATRIP-Del HU/shATRIP-Del HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:`SUV4-20H2`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","SUV4-20H2"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="ATRIP-rescue/EdU-SUV420H2-H3亮度统计（shATRIP rescue）/shATRIP-Del-HU.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-SUV420H2-H3亮度统计（shATRIP rescue）/shATRIP-vector HU/shATRIP-vector HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:`SUV4-20H2`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","SUV4-20H2"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="ATRIP-rescue/EdU-SUV420H2-H3亮度统计（shATRIP rescue）/shATRIP-vector-HU.pdf",height = 2,width = 4)

dat <- read_xlsx("ATRIP-rescue/EdU-SUV420H2-H3亮度统计（shATRIP rescue）/shATRIP-WT HU/shATRIP-WT HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:`SUV4-20H2`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","SUV4-20H2"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="ATRIP-rescue/EdU-SUV420H2-H3亮度统计（shATRIP rescue）/shATRIP-WT-HU.pdf",height = 2,width = 4)

############ Ctrl-HU => EdU-ATRIP ###########
colors19 = c("#b42e20","#a349a4","#313c63")
names(colors19) = c("H3","EdU","ATRIP")
dat <- read_csv("Ctrl-HU/EdU-ATRIP-H3亮度统计 in HeLa WT cell/CTR/EdU-ATRIP CTR.csv") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:ATRIP)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","ATRIP"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="Ctrl-HU/EdU-ATRIP-CTR.pdf",height = 2,width = 4)

dat <- read_csv("Ctrl-HU/EdU-ATRIP-H3亮度统计 in HeLa WT cell/HU/EdU-ATRIP HU.csv") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:ATRIP)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","ATRIP"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="Ctrl-HU/EdU-ATRIP-HU.pdf",height = 2,width = 4)

############ Ctrl-HU => EdU-H4K20me3 ###########
colors19 = c("#b42e20","#a349a4","#5d84a4")
names(colors19) = c("H3","EdU","H4K20me3")
dat <- read_csv("Ctrl-HU/EdU-H4K20me3亮度统计 in HeLa WT cell/CTR/EdU-H4K20me3 CTR.csv") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:H4K20me3)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","H4K20me3"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="Ctrl-HU/EdU-H4K20me3-CTR.pdf",height = 2,width = 4)

setwd("/mnt/sda/Data_ZhuQian/PLA/Figure-LinePlot")
dat <- read_xlsx("Ctrl-HU/EdU-H4K20me3 HU.xlsx") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:H4K20me3)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","H4K20me3"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="Ctrl-HU/EdU-H4K20me3-HU.pdf",height = 2,width = 4)

############ Ctrl-HU => EdU-HP1B ###########
colors19 = c("#b42e20","#a349a4","#7bc7cd")
names(colors19) = c("H3","EdU","HP1B")
dat <- fread("Ctrl-HU/EdU-HP1B-H3亮度统计 in HeLa WT cell/CTR/EdU-HP1B CTR.csv") %>%
  magrittr::set_colnames(c("Distance","H3","EdU","HP1B")) %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:HP1B)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","HP1B"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="Ctrl-HU/EdU-HP1B-CTR.pdf",height = 2,width = 4)

dat <- fread("Ctrl-HU/EdU-HP1B-H3亮度统计 in HeLa WT cell/HU/EdU-HP1B HU.csv") %>%
  magrittr::set_colnames(c("Distance","H3","EdU","HP1B")) %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:HP1B)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","HP1B"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="Ctrl-HU/EdU-HP1B-HU.pdf",height = 2,width = 4)


############ Ctrl-HU => SUV420H2 ###########
colors19 = c("#b42e20","#a349a4","#5d84a4")
names(colors19) = c("H3","EdU","SUV4-20H2")
dat <- fread("Ctrl-HU/EdU-SUV420H2亮度统计 in HeLa WT cell/CTR/EdU-SUV420H2 CTR.csv") %>%
  magrittr::set_colnames(c("Distance","H3","EdU","SUV4-20H2")) %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:`SUV4-20H2`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","SUV4-20H2"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="Ctrl-HU/EdU-SUV4-20H2-CTR.pdf",height = 2,width = 4)

dat <- read_xlsx("Ctrl-HU/EdU-SUV420H2 HU.xlsx") %>%
  magrittr::set_colnames(c("Distance","H3","EdU","SUV4-20H2")) %>%
  magrittr::set_colnames(str_remove_all(colnames(.),"\\(.*") %>% str_remove_all("_.*")) %>%
  gather(group,insensity,H3:`SUV4-20H2`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","SUV4-20H2"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  theme_classic()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    #panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="Ctrl-HU/EdU-SUV4-20H2-HU.pdf",height = 2,width = 4)

############ F2.E ###########
colors19 = c("#b42e20","#995181","#5d84a4")
names(colors19) = c("H3","EdU","SUV4-20H2")
dat <- read_xlsx("F2.E/EdU-SUV420H2-CTR.xlsx") %>% 
  magrittr::set_colnames(c("Distance","EdU","SUV4-20H2")) %>%
  gather(group,insensity,EdU:`SUV4-20H2`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","SUV4-20H2"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="F2.E/EdU-SUV420H2-CTR.pdf",height = 2,width = 4)

dat <- read_xlsx("F2.E/EdU-SUV420H2-HU.xlsx") %>% 
  magrittr::set_colnames(c("Distance","EdU","SUV4-20H2")) %>%
  gather(group,insensity,EdU:`SUV4-20H2`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","SUV4-20H2"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="F2.E/EdU-SUV420H2-HU.pdf",height = 2,width = 4)

############ F2.F ###########
colors19 = c("#b42e20","#995181","#5d84a4")
names(colors19) = c("H3","EdU","H4K20me3")
dat <- read_xlsx("F2.F/EdU-H4K20me3-CTR.xlsx") %>% 
  magrittr::set_colnames(c("Distance","EdU","H4K20me3")) %>%
  gather(group,insensity,EdU:H4K20me3)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","H4K20me3"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="F2.F/EdU-H4K20me3-CTR.pdf",height = 2,width = 4)

dat <- read_xlsx("F2.F/EdU-H4K20me3-HU.xlsx") %>% 
  magrittr::set_colnames(c("Distance","EdU","H4K20me3")) %>%
  gather(group,insensity,EdU:`H4K20me3`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","H4K20me3"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  #geom_line()+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="F2.F/EdU-H4K20me3-HU.pdf",height = 2,width = 4)

############ F2.C ###########
colors19 = c("#b42e20","#995181","#5d84a4")
names(colors19) = c("H3","EdU","ATRIP")
dat <- read_xlsx("F2.C/EdU-ATRIP-CTR.xlsx") %>% 
  magrittr::set_colnames(c("Distance","EdU","ATRIP")) %>%
  gather(group,insensity,EdU:`ATRIP`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","ATRIP"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="F2.C/EdU-ATRIP-CTR.pdf",height = 2,width = 4)

dat <- read_xlsx("F2.C/EdU-ATRIP-HU.xlsx") %>% 
  magrittr::set_colnames(c("Distance","EdU","ATRIP")) %>%
  gather(group,insensity,EdU:`ATRIP`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","ATRIP"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="F2.C/EdU-ATRIP-HU.pdf",height = 2,width = 4)
############ F2.D ###########
colors19 = c("#b42e20","#995181","#5d84a4")
names(colors19) = c("H3","EdU","HP1B")
dat <- read_xlsx("F2.D/EdU-HP1B-CTR.xlsx") %>% 
  magrittr::set_colnames(c("Distance","EdU","HP1B")) %>%
  gather(group,insensity,EdU:`HP1B`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","HP1B"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="F2.D/EdU-HP1B-CTR.pdf",height = 2,width = 4)

dat <- read_xlsx("F2.D/EdU-HP1B-HU.xlsx") %>% 
  magrittr::set_colnames(c("Distance","EdU","HP1B")) %>%
  gather(group,insensity,EdU:`HP1B`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","HP1B"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="F2.D/EdU-HP1B-HU.pdf",height = 2,width = 4)

############ F3.M => H4K20me3 => shATRIP ###########
colors19 = c("#b42e20","#ba5b63",'#377b4c','#7bc7cd','#5d84a4','#4B4B5D')
names(colors19) = c("H3","EdU","ATRIP","HP1B","H4K20me3","SUV420H2")

files = list.files("F3.M/","EdU-H4K20me3.*csv$")
for (file in files) {
  dat <- fread(file.path("F3.M/",file)) %>%
    magrittr::set_colnames(c("Distance","EdU","H4K20me3")) %>%
    gather(group,insensity,EdU:H4K20me3)
  filename = paste(str_replace_all(file," ","-") %>% str_remove_all(".csv"),".pdf",sep="")
  p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","H4K20me3"))),
               aes(Distance,y=insensity,color=group))+
    geom_xspline(spline_shape = 0.6,size=1.2)+
    scale_color_manual("",values = colors19) +
    labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
    ggthemes::theme_few()+
    theme(#legend.position = "bottom",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
      axis.text.y = element_text(size=9),
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size=9),
      panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
      legend.title = element_text(size=9),
      plot.title = element_text(size=9,hjust=0.5),
      legend.text = element_text(size=9))
  ggsave(p1,filename=file.path("F3.M/",filename),height = 2,width = 4)
}

############ F3.M => HP1B => shATRIP ###########
colors19 = c("#b42e20","#ba5b63",'#377b4c','#7bc7cd','#5d84a4','#4B4B5D')
names(colors19) = c("H3","EdU","ATRIP","HP1B","H4K20me3","SUV420H2")

files = list.files("F3.M/","EdU-HP1B.*csv$")
for (file in files) {
  dat <- fread(file.path("F3.M/",file)) %>%
    magrittr::set_colnames(c("Distance","EdU","HP1B")) %>%
    gather(group,insensity,EdU:HP1B)
  filename = paste(str_replace_all(file," ","-") %>% str_remove_all(".csv"),".pdf",sep="")
  p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","HP1B"))),
               aes(Distance,y=insensity,color=group))+
    geom_xspline(spline_shape = 0.6,size=1.2)+
    scale_color_manual("",values = colors19) +
    labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
    ggthemes::theme_few()+
    theme(#legend.position = "bottom",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
      axis.text.y = element_text(size=9),
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size=9),
      panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
      legend.title = element_text(size=9),
      plot.title = element_text(size=9,hjust=0.5),
      legend.text = element_text(size=9))
  ggsave(p1,filename=file.path("F3.M/",filename),height = 2,width = 4)
}

dat <- read_xlsx("F3.M/EdU-HP1B-shATRIP-vector HU.xlsx")
write.csv(dat,file = "F3.M/EdU-HP1B-shATRIP-vector HU.csv",row.names = F,quote = F)

############ F3.M => SUV420H2 => shATRIP ###########
colors19 = c("#b42e20","#ba5b63",'#377b4c','#7bc7cd','#5d84a4','#f5b366')
names(colors19) = c("H3","EdU","ATRIP","HP1B","H4K20me3","SUV420H2")
files = list.files("F3.M/","EdU-SUV420H2.*csv$")
for (file in files) {
  dat <- fread(file.path("F3.M/",file)) %>%
    magrittr::set_colnames(c("Distance","EdU","SUV420H2")) %>%
    gather(group,insensity,EdU:SUV420H2)
  filename = paste(str_replace_all(file," ","-") %>% str_remove_all(".csv"),".pdf",sep="")
  p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","SUV420H2"))),
               aes(Distance,y=insensity,color=group))+
    geom_xspline(spline_shape = 0.6,size=1.2)+
    scale_color_manual("",values = colors19) +
    labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
    ggthemes::theme_few()+
    theme(#legend.position = "bottom",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
      axis.text.y = element_text(size=9),
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size=9),
      panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
      legend.title = element_text(size=9),
      plot.title = element_text(size=9,hjust=0.5),
      legend.text = element_text(size=9))
  ggsave(p1,filename=file.path("F3.M/",filename),height = 2,width = 4)
}

dat <- read_xlsx("F3.M/EdU-SUV420H2-shATRIP-DEL HU.xlsx")
write.csv(dat,file = "F3.M/EdU-SUV420H2-shATRIP-DEL HU.csv",row.names = F,quote = F)

dat <- read_xlsx("F3.M/EdU-SUV420H2-shATRIP-WT HU.xlsx")
write.csv(dat,file = "F3.M/EdU-SUV420H2-shATRIP-WT HU.csv",row.names = F,quote = F)

dat <- read_xlsx("F3.M/EdU-SUV420H2-shATRIP-vector HU.xlsx")
write.csv(dat,file = "F3.M/EdU-SUV420H2-shATRIP-vector HU.csv",row.names = F,quote = F)

############ F3.N => H4K20me3 => shHP1B ###########
colors19 = c("#b42e20","#ba5b63",'#377b4c','#7bc7cd','#5d84a4','#f5b366')
names(colors19) = c("H3","EdU","ATRIP","HP1B","H4K20me3","SUV420H2")
files = list.files("F3.N/","EdU-H4K20me3.*csv$")
for (file in files) {
  dat <- fread(file.path("F3.N/",file)) %>%
    magrittr::set_colnames(c("Distance","EdU","H4K20me3")) %>%
    gather(group,insensity,EdU:H4K20me3)
  filename = paste(str_replace_all(file," ","-") %>% str_remove_all(".csv"),".pdf",sep="")
  p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","H4K20me3"))),
               aes(Distance,y=insensity,color=group))+
    geom_xspline(spline_shape = 0.6,size=1.2)+
    scale_color_manual("",values = colors19) +
    labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
    ggthemes::theme_few()+
    theme(#legend.position = "bottom",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
      axis.text.y = element_text(size=9),
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size=9),
      panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
      legend.title = element_text(size=9),
      plot.title = element_text(size=9,hjust=0.5),
      legend.text = element_text(size=9))
  ggsave(p1,filename=file.path("F3.N/",filename),height = 2,width = 4)
}

############ F3.N => SUV420H2 => shHP1B ###########
colors19 = c("#b42e20","#ba5b63",'#377b4c','#7bc7cd','#5d84a4','#f5b366')
names(colors19) = c("H3","EdU","ATRIP","HP1B","H4K20me3","SUV420H2")
files = list.files("F3.N/","EdU-SUV420H2.*csv$")
for (file in files) {
  dat <- fread(file.path("F3.N/",file)) %>%
    magrittr::set_colnames(c("Distance","EdU","SUV420H2")) %>%
    gather(group,insensity,EdU:SUV420H2)
  filename = paste(str_replace_all(file," ","-") %>% str_remove_all(".csv"),".pdf",sep="")
  p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","SUV420H2"))),
               aes(Distance,y=insensity,color=group))+
    geom_xspline(spline_shape = 0.6,size=1.2)+
    scale_color_manual("",values = colors19) +
    labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
    ggthemes::theme_few()+
    theme(#legend.position = "bottom",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
      axis.text.y = element_text(size=9),
      #axis.ticks.x = element_blank(),
      axis.title = element_text(size=9),
      panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
      legend.title = element_text(size=9),
      plot.title = element_text(size=9,hjust=0.5),
      legend.text = element_text(size=9))
  ggsave(p1,filename=file.path("F3.N/",filename),height = 2,width = 4)
}

################ -------- #####################
setwd("/mnt/sda/Data_ZhuQian/PLA/DNA_Fiber")
############ F5.A #############
dat <- read_xlsx("F5.A/WT siHP1β siSUV4-20H2 Flag-H4K20R.xlsx") %>%
  gather(group,length) %>% mutate(group=if_else(group=="flag-H4K20R","Flag-H4K20R",group))
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.A.L.Wilcox.csv",row.names = F)
p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("WT","siHP1β","siSUV4-20H2","Flag-H4K20R"))),
       aes(group,y=length,color=group))+
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Fork restart\nIdU track length (μm)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,34)+
  annotate("segment", x=1, xend=2, y=26, yend=26,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=26.5,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=29, yend=29,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=29.5,label="****",size=4)+
  annotate("segment", x=1, xend=4, y=32, yend=32,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=32.5,label="****",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
# "WT","siHP1β","siSUV4-20H2","flag-H4K20R"
ggsave(p2,filename="F5.A.L.pdf",height = 2.4,width = 4)

dat <- read_xlsx("F5.A/restart.xlsx") %>% magrittr::set_colnames(c("group","length"))
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.A.R.Wilcox.csv",row.names = F)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("WT","siHP1β","siSUV4-20H2","Flag-H4K20R"))),
       aes(group,y=length,color=group))+
  geom_violin()+
  #geom_boxplot(width=0.2)+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Restart(%)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(25,98)+
  #ggpubr::stat_compare_means(comparisons = list(c("WT","siHP1β"),c("WT","siSUV4-20H2"),c("WT","Flag-H4K20R")),label = "p.signif",method = "wilcox")+
  #geom_signif(comparisons = list(c("WT","siHP1β"),c("WT","siSUV4-20H2"),c("WT","flag-H4K20R")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  annotate("segment", x=1, xend=2, y=88, yend=88,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=88.5,label="***",size=4)+
  annotate("segment", x=1, xend=3, y=92, yend=92,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=92.5,label="****",size=4)+
  annotate("segment", x=1, xend=4, y=96, yend=96,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=96.5,label="***",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
ggsave(p1,filename="F5.A.R.pdf",height = 2.4,width = 3.9)

############ F5.B ###########
dat <- read_xlsx("F5.B/shATRIP rescue.xlsx") %>%
  gather(group,length)
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.B.L.Wilcox.csv",row.names = F)
p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("WT","shATRIP","shATRIP-WT","shATRIP-DEL"))),
             aes(group,y=length,color=group))+
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Fork restart\nIdU track length (μm)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,32)+
  annotate("segment", x=1, xend=2, y=24, yend=24,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=24.5,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=27, yend=27,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=28,label="n.s.",size=3)+
  annotate("segment", x=1, xend=4, y=30, yend=30,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=30.5,label="****",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")

ggsave(p2,filename="F5.B.L.pdf",height = 2.4,width = 4)

dat <- read_xlsx("F5.B/restart.xlsx") %>% magrittr::set_colnames(c("group","length"))
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.B.R.Wilcox.csv",row.names = F)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("WT","shATRIP","shATRIP-WT","shATRIP-DEL"))),
             aes(group,y=length,color=group))+
  geom_violin()+
  #geom_boxplot(width=0.2)+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Restart(%)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(40,95)+
  #ggpubr::stat_compare_means(comparisons = list(c("WT","siHP1β"),c("WT","siSUV4-20H2"),c("WT","Flag-H4K20R")),label = "p.signif",method = "wilcox")+
  #geom_signif(comparisons = list(c("WT","siHP1β"),c("WT","siSUV4-20H2"),c("WT","flag-H4K20R")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  annotate("segment", x=1, xend=2, y=85, yend=85,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=85.5,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=89, yend=89,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=90.5,label="n.s.",size=3)+
  annotate("segment", x=1, xend=4, y=93, yend=93,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=93.5,label="****",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
ggsave(p1,filename="F5.B.R.pdf",height = 2.4,width = 3.9)

############ F5.C ###########
dat <- read_xlsx("F5.C/shHP1β rescue.xlsx") %>%
  gather(group,length)
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.C.L.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(length~group)
write.csv(dat3,file = "F5.C.L.T.csv",row.names = F)

p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("WT","shHP1β","shHP1β-WT","shHP1β-ΔCSD"))),
             aes(group,y=length,color=group))+
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Fork restart\nIdU track length (μm)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,36)+
  annotate("segment", x=1, xend=2, y=28, yend=28,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=28.5,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=31, yend=31,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=32.5,label="n.s.",size=3)+
  annotate("segment", x=1, xend=4, y=34, yend=34,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=34.5,label="****",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
ggsave(p2,filename="F5.C.L.pdf",height = 2.4,width = 4)

dat <- read_xlsx("F5.C/restart.xlsx") %>% magrittr::set_colnames(c("group","length"))
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.C.R.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(length~group)
write.csv(dat3,file = "F5.C.R.T.csv",row.names = F)

p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("WT","shHP1β","shHP1β-WT","shHP1β-ΔCSD"))),
             aes(group,y=length,color=group))+
  geom_violin()+
  #geom_boxplot(width=0.2)+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Restart(%)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(40,99)+
  #ggpubr::stat_compare_means(comparisons = list(c("WT","siHP1β"),c("WT","siSUV4-20H2"),c("WT","Flag-H4K20R")),label = "p.signif",method = "wilcox")+
  #geom_signif(comparisons = list(c("WT","siHP1β"),c("WT","siSUV4-20H2"),c("WT","flag-H4K20R")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  annotate("segment", x=1, xend=2, y=88, yend=88,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=88.5,label="**",size=4)+
  annotate("segment", x=1, xend=3, y=92, yend=92,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=93.5,label="n.s.",size=3)+
  annotate("segment", x=1, xend=4, y=96, yend=96,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=96.5,label="**",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
ggsave(p1,filename="F5.C.R.pdf",height = 2.4,width = 3.9)

############ F5.D ###########
dat <- read_xlsx("F5.D/shSUV4-20H2 rescue.xlsx") %>%
  gather(group,length)
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.D.L.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(length~group)
write.csv(dat3,file = "F5.D.L.T.csv",row.names = F)

p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("WT","shSUV4-20H2","shSUV4-20H2-WT","shSUV4-20H2-Δclamp","shSUV4-20H2-ΔSET"))),
             aes(group,y=length,color=group))+
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Fork restart\nIdU track length (μm)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,42)+
  annotate("segment", x=1, xend=2, y=31, yend=31,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=31.5,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=34, yend=34,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=35.5,label="n.s.",size=3)+
  annotate("segment", x=1, xend=4, y=37, yend=37,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=37.5,label="****",size=4)+
  annotate("segment", x=1, xend=5, y=40, yend=40,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=40.5,label="****",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
ggsave(p2,filename="F5.D.L.pdf",height = 2.4,width = 4.4)

dat <- read_xlsx("F5.D/restart.xlsx") %>% magrittr::set_colnames(c("group","length"))
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.D.R.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(length~group)
write.csv(dat3,file = "F5.D.R.T.csv",row.names = F)

p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("WT","shSUV4-20H2","shSUV4-20H2-WT","shSUV4-20H2-Δclamp","shSUV4-20H2-ΔSET"))),
             aes(group,y=length,color=group))+
  geom_violin()+
  #geom_boxplot(width=0.2)+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Restart(%)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(30,103)+
  #ggpubr::stat_compare_means(comparisons = list(c("WT","siHP1β"),c("WT","siSUV4-20H2"),c("WT","Flag-H4K20R")),label = "p.signif",method = "wilcox")+
  #geom_signif(comparisons = list(c("WT","siHP1β"),c("WT","siSUV4-20H2"),c("WT","flag-H4K20R")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  annotate("segment", x=1, xend=2, y=89, yend=89,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=89.5,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=93, yend=93,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=95,label="n.s.",size=3)+
  annotate("segment", x=1, xend=4, y=97, yend=97,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=97.5,label="***",size=4)+
  annotate("segment", x=1, xend=5, y=101, yend=101,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=101.5,label="**",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
ggsave(p1,filename="F5.D.R.pdf",height = 2.4,width = 4.4)

############ F5.E ###########
dat <- read_xlsx("F5.E/flag-H4K20.xlsx") %>%
  gather(group,length) %>%
  mutate(group=str_replace_all(group,"flag","Flag"))
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.E.L.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(length~group)
write.csv(dat3,file = "F5.E.L.T.csv",row.names = F)

p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("Flag-H4-WT","Flag-H4K20A","Flag-H4K20R","Flag-H4K20Q"))),
             aes(group,y=length,color=group))+
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Fork restart\nIdU track length (μm)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,47)+
  annotate("segment", x=1, xend=2, y=39, yend=39,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=39.5,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=42, yend=42,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=42.5,label="****",size=4)+
  annotate("segment", x=1, xend=4, y=45, yend=45,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=45.5,label="****",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
ggsave(p2,filename="F5.E.L.pdf",height = 2.4,width = 3.9)

dat <- read_xlsx("F5.E/restart.xlsx") %>% magrittr::set_colnames(c("group","length"))
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.E.R.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(length~group)
write.csv(dat3,file = "F5.E.R.T.csv",row.names = F)

p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("Flag-H4","Flag-H4K20A","Flag-H4K20R","Flag-H4K20Q"))),
             aes(group,y=length,color=group))+
  geom_violin()+
  #geom_boxplot(width=0.2)+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Restart(%)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  #ylim(30,103)+
  #ggpubr::stat_compare_means(comparisons = list(c("WT","siHP1β"),c("WT","siSUV4-20H2"),c("WT","Flag-H4K20R")),label = "p.signif",method = "wilcox")+
  #geom_signif(comparisons = list(c("WT","siHP1β"),c("WT","siSUV4-20H2"),c("WT","flag-H4K20R")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  annotate("segment", x=1, xend=2, y=90, yend=90,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=90.5,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=94, yend=94,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=94.5,label="***",size=4)+
  annotate("segment", x=1, xend=4, y=98, yend=98,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=98.5,label="****",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
ggsave(p1,filename="F5.E.R.pdf",height = 2.4,width = 3.9)

############ F5.F ###########
dat <- read_xlsx("F5.F/shATRIP rescue.xlsx") %>%
  gather(group,length) %>%
  mutate(group=str_replace_all(group,"flag","Flag"))
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.F.L.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(length~group)
write.csv(dat3,file = "F5.F.L.T.csv",row.names = F)

p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("WT","shATRIP","shATRIP-WT","shATRIP-DEL"))),
             aes(group,y=length,color=group))+
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Fork restart\nIdU track length (μm)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,39)+
  annotate("segment", x=1, xend=2, y=31, yend=31,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=31.5,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=33, yend=33,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=34,label="n.s.",size=3)+
  annotate("segment", x=1, xend=4, y=37, yend=37,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=37.5,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=31, yend=31,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=31.5,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=31, yend=31,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=31.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=35, yend=35,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=35.5,label="**",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
ggsave(p2,filename="F5.F.pdf",height = 2.4,width = 3.9)

############ F5.G ###########
dat <- read_xlsx("F5.G/shHP1β rescue.xlsx") %>%
  gather(group,length) %>%
  mutate(group=str_replace_all(group,"flag","Flag"))
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.G.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(length~group)
write.csv(dat3,file = "F5.G.T.csv",row.names = F)

p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("WT","shHP1β","shHP1β-WT","shHP1β-ΔCSD"))),
             aes(group,y=length,color=group))+
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Fork restart\nIdU track length (μm)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,34)+
  annotate("segment", x=1, xend=2, y=26, yend=26,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=26.5,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=28, yend=28,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=29,label="n.s.",size=3)+
  annotate("segment", x=1, xend=4, y=32, yend=32,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=32.5,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=26, yend=26,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=26.5,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=26, yend=26,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=26.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=30, yend=30,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=30.5,label="***",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
ggsave(p2,filename="F5.G.pdf",height = 2.4,width = 3.9)

############ F5.H ###########
dat <- read_xlsx("F5.H/shSUV4-20H2 rescue.xlsx") %>%
  gather(group,length) %>%
  mutate(group=str_replace_all(group,"flag","Flag"))
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.H.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(length~group)
write.csv(dat3,file = "F5.H.T.csv",row.names = F)

p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("WT","shSUV4-20H2","shSUV4-20H2-WT","shSUV4-20H2-Δclamp","shSUV4-20H2-ΔSET"))),
             aes(group,y=length,color=group))+
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Fork restart\nIdU track length (μm)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,46)+
  annotate("segment", x=1, xend=2, y=36, yend=36,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=36.5,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=39, yend=39,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=40,label="n.s.",size=3)+
  annotate("segment", x=1, xend=4, y=42, yend=42,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=42.5,label="****",size=4)+
  annotate("segment", x=1, xend=5, y=45, yend=45,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=45.5,label="****",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
ggsave(p2,filename="F5.H.pdf",height = 2.4,width = 4.4)

############ F5.I ###########
dat <- read_xlsx("F5.I/Flag-H4.xlsx") %>%
  gather(group,length) %>%
  mutate(group=str_replace_all(group,"flag","Flag"))
dat3 <- dat %>% rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F5.I.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(length~group)
write.csv(dat3,file = "F5.I.T.csv",row.names = F)

p2 <- ggplot(dat %>% mutate(group=factor(group,levels=c("Flag-H4-WT","Flag-H4K20A","Flag-H4K20R","Flag-H4K20Q"))),
             aes(group,y=length,color=group))+
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Fork restart\nIdU track length (μm)")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,46)+
  annotate("segment", x=1, xend=2, y=39, yend=39,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=39.5,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=42, yend=42,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=42.5,label="****",size=4)+
  annotate("segment", x=1, xend=4, y=45, yend=45,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=45.5,label="****",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")
ggsave(p2,filename="F5.I.pdf",height = 2.4,width = 4)

################ -------- #####################
setwd("/mnt/sda/Data_ZhuQian/PLA/Metaphase")
geom_uperrorbar <- function(mapping = NULL, data = NULL,
                            stat = "identity", position = "identity",
                            ...,
                            na.rm = FALSE,
                            orientation = NA,
                            show.legend = NA,
                            inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomUperrorbar,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      orientation = orientation,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomUperrorbar <- ggproto("GeomUperrorbar", Geom,
                          default_aes = aes(colour = "black", size = 0.5, linetype = 1, width = 0.3,
                                            alpha = NA),
                          
                          draw_key = draw_key_path,
                          
                          required_aes = c("x|y", "ymin|xmin", "ymax|xmax"),
                          
                          setup_params = function(data, params) {
                            GeomLinerange$setup_params(data, params)
                          },
                          
                          extra_params = c("na.rm", "orientation"),
                          
                          setup_data = function(data, params) {
                            data$flipped_aes <- params$flipped_aes
                            data <- flip_data(data, params$flipped_aes)
                            data$width <- data$width %||%
                              params$width %||% (resolution(data$x, FALSE) * 0.9)
                            data <- transform(data,
                                              xmin = x - width / 2, xmax = x + width / 2, width = NULL
                            )
                            flip_data(data, params$flipped_aes)
                          },
                          
                          draw_panel = function(data, panel_params, coord, width = NULL, flipped_aes = FALSE) {
                            data <- flip_data(data, flipped_aes)
                            #x <- as.vector(rbind(data$xmin, data$xmax, NA, data$x,    data$x,    NA, data$xmin, data$xmax))
                            #y <- as.vector(rbind(data$ymax, data$ymax, NA, data$ymax, data$ymin, NA, data$ymin, data$ymin))
                            sel <- data$y < 0 
                            data$ymax[sel] <- data$ymin[sel]
                            x <- as.vector(rbind(data$xmin, data$xmax, NA, data$x,    data$x))
                            y <- as.vector(rbind(data$ymax, data$ymax, NA, data$ymax, data$y))
                            data <- new_data_frame(list(
                              x = x,
                              y = y,
                              colour = rep(data$colour, each = 5),
                              alpha = rep(data$alpha, each = 5),
                              size = rep(data$size, each = 5),
                              linetype = rep(data$linetype, each = 5),
                              group = rep(1:(nrow(data)), each = 5),
                              row.names = 1:(nrow(data) * 5)
                            ))
                            data <- flip_data(data, flipped_aes)
                            GeomPath$draw_panel(data, panel_params, coord)
                          }
)

new_data_frame <- function(x = list(), n = NULL) {
  if (length(x) != 0 && is.null(names(x))) {
    abort("Elements must be named")
  }
  lengths <- vapply(x, length, integer(1))
  if (is.null(n)) {
    n <- if (length(x) == 0 || min(lengths) == 0) 0 else max(lengths)
  }
  for (i in seq_along(x)) {
    if (lengths[i] == n) next
    if (lengths[i] != 1) {
      abort("Elements must equal the number of rows or 1")
    }
    x[[i]] <- rep(x[[i]], n)
  }
  
  class(x) <- "data.frame"
  
  attr(x, "row.names") <- .set_row_names(n)
  x
}
############ F6.C #############
set.seed(123)
dat <- read_xlsx("F6.C WT shATRIP shHP1B shSUV420H2.xlsx") %>% dplyr::select(-...1)%>%
  mutate(Times=paste("R",sample(rep(1:5, each = 20)),sep="")) %>%
  gather(group,length,-Times) %>%
  mutate(Group=cut(length,breaks = c(-Inf,5,15,30,Inf), labels = c("<5","5-15","15-30",">30"), right=FALSE)) %>%
  group_by(group,Group,Times) %>% summarise(Count=n()) %>%
  ungroup() %>% group_by(group,Times) %>% mutate(Total=sum(Count)) %>% ungroup() %>%
  mutate(Proportion=round(Count/Total,10)) %>% ungroup()
write.csv(dat,file = "F6.C.data.csv",row.names = F)

dat2 <- dat %>% group_by(group,Group) %>% summarise(P.m = mean(Proportion),P.sd=sd(Proportion)) %>% ungroup() %>%
  group_by(group) %>%
  mutate(Total=sum(P.m)) %>% ungroup() %>% mutate(Proportion=P.m/Total) %>%
  mutate(Group=factor(Group,levels=rev(c("<5","5-15","15-30",">30")))) %>%
  mutate(group=factor(group,levels=c("WT","shATRIP","shHP1B","shSUV420H2"))) %>%
  arrange(group,Group) %>%
  group_by(group) %>% mutate(xAxis=cumsum(Proportion))
dat2[is.na(dat2)] = 0
write.csv(dat2,file = "F6.C.data-F.csv",row.names = F)

p3 <- dat2 %>% mutate(Group=factor(Group,levels=c("<5","5-15","15-30",">30"))) %>%
  mutate(group=factor(group,levels=c("WT","shATRIP","shHP1B","shSUV420H2"))) %>%
  ggplot(aes(x=group,y=Proportion))+
  geom_bar(aes(fill=Group),stat="identity",width=0.7,size=0.25,alpha=0.7)+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = 'mean_sd', geom = "uperrorbar", colour = "black", width = .3)+
  geom_errorbar(aes(ymin=xAxis, ymax=xAxis+P.sd,color=Group), width=.1)+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="Chromosomal aberrations (%)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c(0,0.2,0.4,0.6,0.8,1))+
  ggtitle("")
ggsave(p3,filename="F6.C.pdf",height = 2.4,width = 3.5)

c("#b1dddd","#81aeac","#598081","#324f4e")

############ F6.D #############
set.seed(123)
dat <- read_xlsx("F6.D Flag-H4 WT K20A K20R K20Q.xlsx") %>% dplyr::select(-...1)%>%
  mutate(Times=paste("R",sample(rep(1:5, each = 20)),sep="")) %>%
  gather(group,length,-Times) %>%
  mutate(Group=cut(length,breaks = c(-Inf,5,15,30,Inf), labels = c("<5","5-15","15-30",">30"), right=FALSE)) %>%
  group_by(group,Group,Times) %>% summarise(Count=n()) %>%
  ungroup() %>% group_by(group,Times) %>% mutate(Total=sum(Count)) %>% ungroup() %>%
  mutate(Proportion=round(Count/Total,10)) %>% ungroup()
write.csv(dat,file = "F6.D.data.csv",row.names = F)

dat2 <- dat %>% group_by(group,Group) %>% summarise(P.m = mean(Proportion),P.sd=sd(Proportion)) %>% ungroup() %>%
  group_by(group) %>%
  mutate(Total=sum(P.m)) %>% ungroup() %>% mutate(Proportion=P.m/Total) %>%
  mutate(Group=factor(Group,levels=rev(c("<5","5-15","15-30",">30")))) %>%
  mutate(group=factor(group,levels=c("Flag-H4 WT","Flag-H4K20A","Flag-H4K20Q","Flag-H4K20R"))) %>%
  arrange(group,Group) %>%
  group_by(group) %>% mutate(xAxis=cumsum(Proportion))
dat2[is.na(dat2)] = 0
write.csv(dat2,file = "F6.D.data-F.csv",row.names = F)

p3 <- dat2 %>% mutate(Group=factor(Group,levels=c("<5","5-15","15-30",">30"))) %>%
  mutate(group=factor(group,levels=c("Flag-H4 WT","Flag-H4K20A","Flag-H4K20Q","Flag-H4K20R"))) %>%
  ggplot(aes(x=group,y=Proportion))+
  geom_bar(aes(fill=Group),stat="identity",width=0.7,size=0.25,alpha=0.7)+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = 'mean_sd', geom = "uperrorbar", colour = "black", width = .3)+
  geom_errorbar(aes(ymin=xAxis, ymax=xAxis+P.sd,color=Group), width=.1)+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="Chromosomal aberrations (%)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c(0,0.2,0.4,0.6,0.8,1))+
  ggtitle("")
ggsave(p3,filename="F6.D.pdf",height = 2.4,width = 3.5)

############ F6.E #############
set.seed(123)
dat <- read_xlsx("F6.E shATRIP rescue.xlsx") %>% dplyr::select(-...1)%>%
  mutate(Times=paste("R",sample(rep(1:5, each = 20)),sep="")) %>%
  gather(group,length,-Times) %>%
  mutate(Group=cut(length,breaks = c(-Inf,5,15,30,Inf), labels = c("<5","5-15","15-30",">30"), right=FALSE)) %>%
  group_by(group,Group,Times) %>% summarise(Count=n()) %>%
  ungroup() %>% group_by(group,Times) %>% mutate(Total=sum(Count)) %>% ungroup() %>%
  mutate(Proportion=round(Count/Total,10)) %>% ungroup()
write.csv(dat,file = "F6.E.data.csv",row.names = F)

dat2 <- dat %>% group_by(group,Group) %>% summarise(P.m = mean(Proportion),P.sd=sd(Proportion)) %>% ungroup() %>%
  group_by(group) %>%
  mutate(Total=sum(P.m)) %>% ungroup() %>% mutate(Proportion=P.m/Total) %>%
  mutate(Group=factor(Group,levels=rev(c("<5","5-15","15-30",">30")))) %>%
  mutate(group=factor(group,levels=c("WT","shATRIP-Vector","shATRIP-WT","shATRIP-Del"))) %>%
  arrange(group,Group) %>%
  group_by(group) %>% mutate(xAxis=cumsum(Proportion))
dat2[is.na(dat2)] = 0
write.csv(dat2,file = "F6.E.data-F.csv",row.names = F)

p3 <- dat2 %>% mutate(Group=factor(Group,levels=c("<5","5-15","15-30",">30"))) %>%
  mutate(group=factor(group,levels=c("WT","shATRIP-Vector","shATRIP-WT","shATRIP-Del"))) %>%
  ggplot(aes(x=group,y=Proportion))+
  geom_bar(aes(fill=Group),stat="identity",width=0.7,size=0.25,alpha=0.7)+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = 'mean_sd', geom = "uperrorbar", colour = "black", width = .3)+
  geom_errorbar(aes(ymin=xAxis, ymax=xAxis+P.sd,color=Group), width=.1)+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="Chromosomal aberrations (%)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1[2:4])+
  scale_color_manual("",values = colorlist1[2:4])+
  guides(color = guide_legend(override.aes = list(size=3)))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c(0,0.2,0.4,0.6,0.8,1))+
  ggtitle("")
ggsave(p3,filename="F6.E.pdf",height = 2.4,width = 3.5)
############ F6.F #############
set.seed(123)
dat <- read_xlsx("F6.F shHP1B rescue.xlsx") %>% dplyr::select(-...1)%>%
  mutate(Times=paste("R",sample(rep(1:5, each = 20)),sep="")) %>%
  gather(group,length,-Times) %>%
  mutate(Group=cut(length,breaks = c(-Inf,5,15,30,Inf), labels = c("<5","5-15","15-30",">30"), right=FALSE)) %>%
  group_by(group,Group,Times) %>% summarise(Count=n()) %>%
  ungroup() %>% group_by(group,Times) %>% mutate(Total=sum(Count)) %>% ungroup() %>%
  mutate(Proportion=round(Count/Total,10)) %>% ungroup()
write.csv(dat,file = "F6.F.data.csv",row.names = F)

dat2 <- dat %>% group_by(group,Group) %>% summarise(P.m = mean(Proportion),P.sd=sd(Proportion)) %>% ungroup() %>%
  group_by(group) %>%
  mutate(Total=sum(P.m)) %>% ungroup() %>% mutate(Proportion=P.m/Total) %>%
  mutate(Group=factor(Group,levels=rev(c("<5","5-15","15-30",">30")))) %>%
  mutate(group=factor(group,levels=c("WT","shHP1B-Vector","shHP1B-WT","shHP1B-ΔCSD"))) %>%
  arrange(group,Group) %>%
  group_by(group) %>% mutate(xAxis=cumsum(Proportion))
dat2[is.na(dat2)] = 0
write.csv(dat2,file = "F6.F.data-F.csv",row.names = F)

p3 <- dat2 %>% mutate(Group=factor(Group,levels=c("<5","5-15","15-30",">30"))) %>%
  mutate(group=factor(group,levels=c("WT","shHP1B-Vector","shHP1B-WT","shHP1B-ΔCSD"))) %>%
  ggplot(aes(x=group,y=Proportion))+
  geom_bar(aes(fill=Group),stat="identity",width=0.7,size=0.25,alpha=0.7)+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = 'mean_sd', geom = "uperrorbar", colour = "black", width = .3)+
  geom_errorbar(aes(ymin=xAxis, ymax=xAxis+P.sd,color=Group), width=.1)+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="Chromosomal aberrations (%)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1[2:4])+
  scale_color_manual("",values = colorlist1[2:4])+
  guides(color = guide_legend(override.aes = list(size=3)))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c(0,0.2,0.4,0.6,0.8,1))+
  ggtitle("")
ggsave(p3,filename="F6.F.pdf",height = 2.4,width = 3.5)

############ F6.G #############
set.seed(123)
dat <- read_xlsx("F6.G shSUV420H2 rescue.xlsx") %>% dplyr::select(-...1)%>%
  mutate(Times=paste("R",sample(rep(1:5, each = 20)),sep="")) %>%
  gather(group,length,-Times) %>%
  mutate(Group=cut(length,breaks = c(-Inf,5,15,30,Inf), labels = c("<5","5-15","15-30",">30"), right=FALSE)) %>%
  dplyr::group_by(group,Group,Times) %>% dplyr::summarise(Count=n()) %>%
  dplyr::ungroup() %>% dplyr::group_by(group,Times) %>% dplyr::mutate(Total=sum(Count)) %>% dplyr::ungroup() %>%
  dplyr::mutate(Proportion=round(Count/Total,10)) %>% dplyr::ungroup()
write.csv(dat,file = "F6.G.data.csv",row.names = F)

dat2 <- dat %>% dplyr::group_by(group,Group) %>% dplyr::summarise(P.m = mean(Proportion),P.sd=sd(Proportion)) %>% dplyr::ungroup() %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(Total=sum(P.m)) %>% dplyr::ungroup() %>% dplyr::mutate(Proportion=P.m/Total) %>%
  dplyr::mutate(Group=factor(Group,levels=rev(c("<5","5-15","15-30",">30")))) %>%
  dplyr::mutate(group=factor(group,levels=c("WT","shSUV420H2-vector","shSUV420H2-WT","shSUV420H2-Δclamp","shSUV420H2-ΔSET"))) %>%
  dplyr::arrange(group,Group) %>%
  dplyr::group_by(group) %>%dplyr:: mutate(xAxis=cumsum(Proportion))
dat2[is.na(dat2)] = 0
write.csv(dat2,file = "F6.G.data-F.csv",row.names = F)

p3 <- dat2 %>% mutate(Group=factor(Group,levels=c("<5","5-15","15-30",">30"))) %>%
  mutate(group=factor(group,levels=c("WT","shSUV420H2-vector","shSUV420H2-WT","shSUV420H2-Δclamp","shSUV420H2-ΔSET"))) %>%
  ggplot(aes(x=group,y=Proportion))+
  geom_bar(aes(fill=Group),stat="identity",width=0.7,size=0.25,alpha=0.7)+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = 'mean_sd', geom = "uperrorbar", colour = "black", width = .3)+
  geom_errorbar(aes(ymin=xAxis, ymax=xAxis+P.sd,color=Group), width=.1)+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="Chromosomal aberrations (%)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c(0,0.2,0.4,0.6,0.8,1))+
  ggtitle("")
ggsave(p3,filename="F6.G.pdf",height = 2.4,width = 3.5)

############ F6.A ##############
setwd("/mnt/sda/Data_ZhuQian/PLA/Micronucleus")
dat <- read_xlsx("F6.A WT shATRIP shHP1B shSUV420H2 Micronucleus per cell - 副本.xlsx")
colnames(dat)[1] = "group"
dat2 <- dat %>% gather(Group,length,-group) %>% 
  mutate(group=factor(group,levels=c("WT","shATRIP","shHP1β","shSUV420H2")))

dat3 <- dat2 %>% group_by(Group) %>%
  rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F6.A.W.csv",row.names = F)

dat4 <- dat2 %>% group_by(Group) %>%
  rstatix::t_test(length~group)
write.csv(dat4,file = "F6.A.T.csv",row.names = F)

p5 <- dat2 %>% 
  mutate(group=factor(group,levels=c("WT","shATRIP","shHP1β","shSUV420H2"))) %>%
  ggplot(aes(x=Group,y=length,fill=group))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.3,jitter.height = 0.001),size=2)+
  stat_summary(fun.data = 'mean_se', geom = "uperrorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="Micronucleus per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")+
  annotate("segment", x=0.7, xend=0.9, y=0.6, yend=0.6,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=0.8,y=0.62,label="****",size=4)+
  annotate("segment", x=0.7, xend=1.1, y=0.68, yend=0.68,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=0.9,y=0.7,label="****",size=4)+
  annotate("segment", x=0.7, xend=1.3, y=0.76, yend=0.76,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1,y=0.78,label="**",size=4)+
  annotate("segment", x=1.7, xend=1.9, y=0.6, yend=0.6,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.8,y=0.62,label="***",size=4)+
  annotate("segment", x=1.7, xend=2.1, y=0.68, yend=0.68,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.9,y=0.7,label="****",size=4)+
  annotate("segment", x=1.7, xend=2.3, y=0.76, yend=0.76,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=0.78,label="***",size=4)
ggsave(p5,filename="F6.A.pdf",height = 2.4,width=4)

############ F6.B ##############
setwd("/mnt/sda/Data_ZhuQian/PLA/Micronucleus")

dat <- read_xlsx("F6.B Flag-H4  Micronucleus per cell.xlsx")
colnames(dat)[1] = "group"
dat2 <- dat %>% gather(Group,length,-group) %>% 
  mutate(group=factor(group,levels=c("Flag-H4","Flag-H4K20A","Flag-H4K20R","Flag-H4K20Q")))

dat3 <- dat2 %>% group_by(Group) %>%
  rstatix::wilcox_test(length~group)
write.csv(dat3,file = "F6.B.W.csv",row.names = F)

dat4 <- dat2 %>% group_by(Group) %>%
  rstatix::t_test(length~group)
write.csv(dat4,file = "F6.B.T.csv",row.names = F)

p5 <- dat2 %>% 
  mutate(group=factor(group,levels=c("Flag-H4","Flag-H4K20A","Flag-H4K20R","Flag-H4K20Q"))) %>%
  ggplot(aes(x=Group,y=length,fill=group))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.3,jitter.height = 0.000001),size=2)+
  stat_summary(fun.data = 'mean_se', geom = "uperrorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="Micronucleus per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")+
  annotate("segment", x=0.7, xend=0.9, y=0.6, yend=0.6,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=0.8,y=0.62,label="**",size=4)+
  annotate("segment", x=0.7, xend=1.1, y=0.68, yend=0.68,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=0.9,y=0.7,label="***",size=4)+
  annotate("segment", x=0.7, xend=1.3, y=0.76, yend=0.76,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1,y=0.78,label="***",size=4)+
  annotate("segment", x=1.7, xend=1.9, y=0.6, yend=0.6,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.8,y=0.62,label="****",size=4)+
  annotate("segment", x=1.7, xend=2.1, y=0.68, yend=0.68,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.9,y=0.7,label="**",size=4)+
  annotate("segment", x=1.7, xend=2.3, y=0.76, yend=0.76,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=0.78,label="**",size=4)+
  ylim(0,0.82)
ggsave(p5,filename="F6.B.pdf",height = 2.4,width=4)

################ -- Clone -- ####################
setwd("/mnt/sda/Data_ZhuQian/PLA/Clones/")
############ F7.C ############
dat2 <- read_xlsx("siA_B_K-HU.xlsx") %>% magrittr::set_colnames(c("NC","siATRIP","siHP1","siKMT5C","H4K20R")) %>%
  gather(group,length)
dat2 %>% rstatix::t_test(length~group)
p51 <- dat2 %>% 
  mutate(group=factor(group,levels=c("NC","siATRIP","siHP1","siKMT5C","H4K20R"))) %>%
  ggplot(aes(x=group,y=length,fill=group))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.3,jitter.height = 0.001),size=2)+
  stat_summary(fun.data = 'mean_se', geom = "uperrorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="No. of clones")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")+
  ylim(0,230)
ggsave(p51,filename="F7.C.pdf",height = 2.4,width=4)

############ F7.D ############
dat2 <- read_xlsx("siArescue_HU.xlsx") %>% magrittr::set_colnames(c("NC","NC_HU","Vec_HU","WT_HU","Del_HU")) %>%
  gather(group,length)
dat2 %>% rstatix::t_test(length~group)
p521 <- dat2 %>% 
  mutate(group=factor(group,levels=c("NC","NC_HU","Vec_HU","WT_HU","Del_HU"))) %>%
  ggplot(aes(x=group,y=length,fill=group))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.3,jitter.height = 0.001),size=2)+
  stat_summary(fun.data = 'mean_se', geom = "uperrorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="No. of clones")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")+
  ylim(0,260)
ggsave(p521,filename="F7.D.pdf",height = 2.4,width=3.5)

############ F7.h ############
dat2 <- read_xlsx("ctr_siATRIP_siHP1B-siKMT5C+Pi.xlsx") %>% magrittr::set_colnames(c("NC","siATRIP","siHP1","siKMT5C")) %>%
  gather(group,length)
dat2 %>% rstatix::t_test(length~group)
p5211 <- dat2 %>% 
  mutate(group=factor(group,levels=c("NC","siATRIP","siHP1","siKMT5C"))) %>%
  ggplot(aes(x=group,y=length,fill=group))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.3,jitter.height = 0.001),size=2)+
  stat_summary(fun.data = 'mean_se', geom = "uperrorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="No. of clones")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")+
  ylim(0,220)
ggsave(p5211,filename="F7.h.pdf",height = 2.4,width=4)

############ F7.i ############
dat2 <- read_xlsx("CTR_A196_Pi_A196+Pi.xlsx") %>% magrittr::set_colnames(c("NC","A196","PARPi","A196_PARPi")) %>%
  gather(group,length)
dat2 %>% rstatix::t_test(length~group)
p5211 <- dat2 %>% 
  mutate(group=factor(group,levels=c("NC","A196","PARPi","A196_PARPi"))) %>%
  ggplot(aes(x=group,y=length,fill=group))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.3,jitter.height = 0.001),size=2)+
  stat_summary(fun.data = 'mean_se', geom = "uperrorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="No. of clones")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")+
  ylim(0,250)
ggsave(p5211,filename="F7.i.pdf",height = 2.4,width=4)

############ F7.j ############
dat2 <- read_xlsx("ctr_Pi_siA-VEC_siA-WT_siA-DEL.xlsx") %>% magrittr::set_colnames(c("NC","NC_Pi","Vec_Pi","WT_Pi","Del_Pi")) %>%
  gather(group,length)
dat2 %>% rstatix::t_test(length~group)
p5211 <- dat2 %>% 
  mutate(group=factor(group,levels=c("NC","NC_Pi","Vec_Pi","WT_Pi","Del_Pi"))) %>%
  ggplot(aes(x=group,y=length,fill=group))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.3,jitter.height = 0.001),size=2)+
  stat_summary(fun.data = 'mean_se', geom = "uperrorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="No. of clones")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")+
  ylim(0,240)
ggsave(p5211,filename="F7.j.pdf",height = 2.4,width=4)


################ -- Clone -2 -- ####################
setwd("/mnt/sda/Data_ZhuQian/PLA/Clones/Clones-2/")
############ F7.C ############
dat2 <- read_xlsx("siA_B_K-HU.xlsx") %>% magrittr::set_colnames(c("NC","siATRIP","siHP1","siKMT5C","H4K20R")) %>%
  gather(group,length)
dat2 %>% rstatix::t_test(length~group)
p51 <- dat2 %>% 
  mutate(group=factor(group,levels=c("NC","siATRIP","siHP1","siKMT5C","H4K20R"))) %>%
  ggplot(aes(x=group,y=length,fill=group))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.3,jitter.height = 0.001),size=2)+
  stat_summary(fun.data = 'mean_se', geom = "uperrorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="No. of clones")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")+
  ylim(0,400)
ggsave(p51,filename="F7.C.pdf",height = 2.4,width=4)

############ F7.D ############
dat2 <- read_xlsx("siArescue_HU.xlsx") %>% magrittr::set_colnames(c("NC","NC_HU","Vec_HU","WT_HU","Del_HU")) %>%
  gather(group,length)
dat2 %>% rstatix::t_test(length~group)
p521 <- dat2 %>% 
  mutate(group=factor(group,levels=c("NC","NC_HU","Vec_HU","WT_HU","Del_HU"))) %>%
  ggplot(aes(x=group,y=length,fill=group))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.3,jitter.height = 0.001),size=2)+
  stat_summary(fun.data = 'mean_se', geom = "uperrorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="No. of clones")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")+
  ylim(0,500)
ggsave(p521,filename="F7.D.pdf",height = 2.4,width=3.5)

############ F7.h ############
dat2 <- read_xlsx("ctr_siATRIP_siHP1B-siKMT5C+Pi.xlsx") %>% magrittr::set_colnames(c("NC","siATRIP","siHP1","siKMT5C")) %>%
  gather(group,length)
dat2 %>% rstatix::t_test(length~group)
p5211 <- dat2 %>% 
  mutate(group=factor(group,levels=c("NC","siATRIP","siHP1","siKMT5C"))) %>%
  ggplot(aes(x=group,y=length,fill=group))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.3,jitter.height = 0.001),size=2)+
  stat_summary(fun.data = 'mean_se', geom = "uperrorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="No. of clones")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")+
  ylim(0,400)
ggsave(p5211,filename="F7.h.pdf",height = 2.4,width=4)

############ F7.i ############
dat2 <- read_xlsx("CTR_A196_Pi_A196+Pi.xlsx") %>% magrittr::set_colnames(c("NC","A196","PARPi","A196_PARPi")) %>%
  gather(group,length)
dat2 %>% rstatix::t_test(length~group)
p5211 <- dat2 %>% 
  mutate(group=factor(group,levels=c("NC","A196","PARPi","A196_PARPi"))) %>%
  ggplot(aes(x=group,y=length,fill=group))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.3,jitter.height = 0.001),size=2)+
  stat_summary(fun.data = 'mean_se', geom = "uperrorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="No. of clones")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")+
  ylim(0,500)
ggsave(p5211,filename="F7.i.pdf",height = 2.4,width=4)

############ F7.j ############
dat2 <- read_xlsx("ctr_Pi_siA-VEC_siA-WT_siA-DEL.xlsx") %>% magrittr::set_colnames(c("NC","NC_Pi","Vec_Pi","WT_Pi","Del_Pi")) %>%
  gather(group,length)
dat2 %>% rstatix::t_test(length~group)
p5211 <- dat2 %>% 
  mutate(group=factor(group,levels=c("NC","NC_Pi","Vec_Pi","WT_Pi","Del_Pi"))) %>%
  ggplot(aes(x=group,y=length,fill=group))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.3,jitter.height = 0.001),size=2)+
  stat_summary(fun.data = 'mean_se', geom = "uperrorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  #scale_y_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  labs(x=NULL,y="No. of clones")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks = element_text(color="black"),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("")+
  ylim(0,500)
ggsave(p5211,filename="F7.j.pdf",height = 2.4,width=4)



############### Figure => Bioinformatics  ###########
############ Figure 1 => Ctrl HU ###########
######## MNase ########
######## ATAC ########
######## MINCE #######
# /mnt/sda/Data_ZhuQian/MINCE-seq-20240905-Second/IP-Del
setwd("/mnt/sda/Data_ZhuQian/MINCE-seq-20240905-Second/IP-Del")
dat1 <- fread("GenomeWindow-2k-CtrvsHU.DESeq2.Up.txt.bed")
dat2 <- fread("GenomeWindow-2k-CtrvsHU.DESeq2.Down.txt.bed")
# HU vs Ctrl
dat <- data.frame(Count=c(nrow(dat1),nrow(dat2)),Group=c("Decreased","Increased"))
write.csv(dat,"GenomeWindow-2k-CtrvsHU.DESeq2-Stat.csv",row.names = F)
colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c','#7bc7cd','#5d84a4','#4B4B5D',"#EC7232")
colors <- c(
  "Decreased" = '#313c63',
  "Increased" = '#b42e20'
)

library(ggforce)
p <- ggplot(dat) +
  geom_link(aes(x = 0, y = Group,
                xend = Count, yend = Group,
                alpha = after_stat(index),
                color = Group,
                size = after_stat(index)),
            n = 500, show.legend = T)
p

p1 <- p +
  geom_point(aes(x = Count,y = Group), color = "black", fill = "white",size = 3,shape = 21) +
  geom_text(aes(x = Count, y = Group,label=Count), size=3, nudge_x=0.05) +
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        legend.text = element_text(size=9)) +
  xlab("No. of 2k windows") + ylab("") +
  scale_color_manual(values = colors)+
  xlim(0,5000)
p1
ggsave(p1,filename = "GenomeWindow-2k-CtrvsHU.DESeq2-Count.pdf", width = 2, height = 1.5)

#### 5k
dat <- fread("IP-CtrHU.DE.Windows/CtrlHU-IZ-SLOP-2Kwindows-COUNT.CSV")
# HU vs Ctrl
dat <- dat %>% filter(Xaxis=="5k") %>% mutate(Group=if_else(Group=="Up","Increased","Decreased"))
#write.csv(dat,"GenomeWindow-2k-CtrvsHU.DESeq2-Stat.csv",row.names = F)
colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c','#7bc7cd','#5d84a4','#4B4B5D',"#EC7232")
colors <- c(
  "Decreased" = '#313c63',
  "Increased" = '#b42e20'
)

library(ggforce)
p <- ggplot(dat) +
  geom_link(aes(x = 0, y = Group,
                xend = Count, yend = Group,
                alpha = after_stat(index),
                color = Group,
                size = after_stat(index)),
            n = 500, show.legend = T)
p

p1 <- p +
  geom_point(aes(x = Count,y = Group), color = "black", fill = "white",size = 3,shape = 21) +
  geom_text(aes(x = Count, y = Group,label=Count), size=3, nudge_x=0.05) +
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        legend.text = element_text(size=9)) +
  xlab("No. of 2k windows") + ylab("") +
  scale_color_manual(values = colors)+
  xlim(0,300)
p1
ggsave(p1,filename = "GenomeWindow-2k-CtrvsHU.DESeq2-Count-5K-Slop.pdf", width = 2, height = 1.5)

### 50 k
dat <- fread("IP-CtrHU.DE.Windows/CtrlHU-IZ-SLOP-2Kwindows-COUNT.CSV")
# HU vs Ctrl
dat <- dat %>% filter(Xaxis=="50k") %>% mutate(Group=if_else(Group=="Up","Increased","Decreased"))
#write.csv(dat,"GenomeWindow-2k-CtrvsHU.DESeq2-Stat.csv",row.names = F)
colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c','#7bc7cd','#5d84a4','#4B4B5D',"#EC7232")
colors <- c(
  "Decreased" = '#313c63',
  "Increased" = '#b42e20'
)

library(ggforce)
p <- ggplot(dat) +
  geom_link(aes(x = 0, y = Group,
                xend = Count, yend = Group,
                alpha = after_stat(index),
                color = Group,
                size = after_stat(index)),
            n = 500, show.legend = T)
p

p1 <- p +
  geom_point(aes(x = Count,y = Group), color = "black", fill = "white",size = 3,shape = 21) +
  geom_text(aes(x = Count, y = Group,label=Count), size=3, nudge_x=0.05) +
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        legend.text = element_text(size=9)) +
  xlab("No. of 2k windows") + ylab("") +
  scale_color_manual(values = colors)+
  xlim(0,1100)
p1
ggsave(p1,filename = "GenomeWindow-2k-CtrvsHU.DESeq2-Count-50K-Slop.pdf", width = 2, height = 1.5)

######## repliATAC #########
setwd("/mnt/sda/Data_ZhuQian/Repli-ATAC-20240905/IP-ClassI/macs3")
##
#/mnt/sda/Data_ZhuQian/Repli-ATAC-20240905/IP-ClassI/macs3

######## MINCE =》 Signal #######
setwd("/mnt/sda/Data_ZhuQian/MINCE-seq-20240905-Second/IP-Del/IP-CtrHU.DE.Windows")
dat <- fread("2Kwindows.Up-IZ.50000.bed-CtrlHU-2kWindows.SR_matrix.gz")
dat <- dat[,-c(1:6)]
dat2 <- colMeans(dat)
dat3 <- data.frame(RPKM=dat2[c(101:200,401:500)],group=rep(c("Ctrl","HU"),c(100,100)))

p1 <- dat3 %>% ggplot(aes(x=group,y=RPKM))+
  geom_violin(aes(fill=group),alpha=0.5)+
  geom_quasirandom(aes(color=group),size = 1, groupOnX = FALSE,alpha=0.4) + # 使用更紧凑排列的蜂群点
  geom_boxplot(width=0.2,color="grey90")+
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "RPKM")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("IZ 2k windows 50kb")+
  annotate("segment", x=1, xend=2, y=0.3, yend=0.3,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=0.31,label="****",size=4)+
  ylim(0,0.34)
ggsave(p1,filename="ggplot2-2Kwindows.Up-IZ.50000.bed-CtrlHU-2kWindows.SR.pdf",height = 2,width = 1.5)

######## MINCE =》 Signal =》 IZ complement ##########
#### 10000000 ####
setwd("/mnt/sda/Data_ZhuQian/MINCE-seq-20240905-Second/IP-Del")
dat <- fread("complement-IZ.10000000.bed.merge.center.matrix.gz")

dat <- dat[,-c(1:6)]
dat2 <- colMeans(dat)
dat3 <- data.frame(RPKM=dat2,group=rep(c("Ctrl","HU"),c(100,100)))
colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c','#7bc7cd','#5d84a4','#4B4B5D',"#EC7232")
p1 <- dat3 %>% ggplot(aes(x=group,y=RPKM))+
  geom_violin(aes(fill=group),alpha=0.5)+
  geom_quasirandom(aes(color=group),size = 1, groupOnX = FALSE,alpha=0.4) + # 使用更紧凑排列的蜂群点
  geom_boxplot(width=0.2,color="grey90")+
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "RPKM")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("complement-IZ.10000000")+
  #annotate("segment", x=1, xend=2, y=0.3, yend=0.3,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  #annotate("text",x=1.5,y=0.31,label="****",size=4)+
  ggpubr::stat_compare_means()
ggsave(p1,filename="complement-IZ.10000000.bed.merge.center.SR.pdf",height = 2,width = 1.5)

#### 1000000 ####
setwd("/mnt/sda/Data_ZhuQian/MINCE-seq-20240905-Second/IP-Del")
dat <- fread("complement-IZ.1000000.bed.merge.center.matrix.gz")

dat <- dat[,-c(1:6)]
dat2 <- colMeans(dat)
dat3 <- data.frame(RPKM=dat2,group=rep(c("Ctrl","HU"),c(100,100)))
colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c','#7bc7cd','#5d84a4','#4B4B5D',"#EC7232")
p1 <- dat3 %>% ggplot(aes(x=group,y=RPKM))+
  geom_violin(aes(fill=group),alpha=0.5)+
  geom_quasirandom(aes(color=group),size = 1, groupOnX = FALSE,alpha=0.4) + # 使用更紧凑排列的蜂群点
  geom_boxplot(width=0.2,color="grey90")+
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "RPKM")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("complement-IZ.1000000")+
  #annotate("segment", x=1, xend=2, y=0.3, yend=0.3,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  #annotate("text",x=1.5,y=0.31,label="****",size=4)+
  ggpubr::stat_compare_means()
ggsave(p1,filename="complement-IZ.1000000.bed.merge.center.SR.pdf",height = 2,width = 1.5)
#### 100000000 ####
setwd("/mnt/sda/Data_ZhuQian/MINCE-seq-20240905-Second/IP-Del")
dat <- fread("complement-IZ.100000000.bed.merge.center.matrix.gz")

dat <- dat[,-c(1:6)]
dat2 <- colMeans(dat)
dat3 <- data.frame(RPKM=dat2,group=rep(c("Ctrl","HU"),c(100,100)))
colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c','#7bc7cd','#5d84a4','#4B4B5D',"#EC7232")
p1 <- dat3 %>% ggplot(aes(x=group,y=RPKM))+
  geom_violin(aes(fill=group),alpha=0.5)+
  geom_quasirandom(aes(color=group),size = 1, groupOnX = FALSE,alpha=0.4) + # 使用更紧凑排列的蜂群点
  geom_boxplot(width=0.2,color="grey90")+
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "RPKM")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("complement-IZ.100000000")+
  #annotate("segment", x=1, xend=2, y=0.3, yend=0.3,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  #annotate("text",x=1.5,y=0.31,label="****",size=4)+
  ggpubr::stat_compare_means()
ggsave(p1,filename="complement-IZ.100000000.bed.merge.center.SR.pdf",height = 2,width = 1.5)
######## repliATAC =》 Signal #########
setwd("/mnt/sda/Data_ZhuQian/Repli-ATAC-20240905/IP-ClassI/macs3")
dat <- fread("IZ50K-CD-IP-CtrHUC.merge.matrix.gz")
dat <- dat[,-c(1:6)]
dat2 <- colMeans(dat)
dat3 <- data.frame(RPKM=dat2[c(1:200,201:400)],group=rep(c("Ctrl","HU"),c(200,200)))

p1 <- dat3 %>% filter(RPKM > 0) %>%
  ggplot(aes(x=group,y=RPKM))+
  geom_violin(aes(fill=group),alpha=0.5)+
  geom_quasirandom(aes(color=group),size = 1, groupOnX = FALSE,alpha=0.4) + # 使用更紧凑排列的蜂群点
  geom_boxplot(width=0.05,color="grey90")+
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "RPKM")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggpubr::stat_compare_means()+
  ggtitle("IZ 2k windows 50kb")+
  annotate("segment", x=1, xend=2, y=0.3, yend=0.3,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=0.31,label="****",size=4)+
  #ylim(0,0.34)
ggsave(p1,filename="ggplot2-2Kwindows.Up-IZ.50000.bed-CtrlHU-2kWindows.SR.pdf",height = 2,width = 1.5)

######## repliATAC =》 Signal =》 IZ complement ##########
#### 1000000
setwd("/mnt/sda/Data_ZhuQian/Repli-ATAC-20240905/IP-ClassI/macs3")
dat <- fread("complement-IZ.1000000.bed.merge.center.matrix.gz")

dat <- dat[,-c(1:6)]
dat[is.na(dat)]=0
dat2 <- colMeans(dat)
dat3 <- data.frame(RPKM=dat2,group=rep(c("Ctrl","HU"),c(100,100))) %>% filter(RPKM <= 0.5)

colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c','#7bc7cd','#5d84a4','#4B4B5D',"#EC7232")
p1 <- dat3 %>% ggplot(aes(x=group,y=RPKM))+
  geom_violin(aes(fill=group),alpha=0.5)+
  geom_quasirandom(aes(color=group),size = 1, groupOnX = FALSE,alpha=0.4) + # 使用更紧凑排列的蜂群点
  geom_boxplot(width=0.1,color="grey90")+
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "RPKM")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("complement-IZ.1000000")+
  #annotate("segment", x=1, xend=2, y=0.3, yend=0.3,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  #annotate("text",x=1.5,y=0.31,label="****",size=4)+
  ggpubr::stat_compare_means()
ggsave(p1,filename="complement-IZ.1000000.bed.merge.center.SR.pdf",height = 2,width = 1.5)

############ Figure 3 ###########
######## ABS - MINCE-seq => Sinal ########
setwd("/mnt/sda/Data_ZhuQian/mince-seq-20240917-Third")
dat <- fread("MINCE-GenomeWindow-2k-AvsHU.DESeq2.Up.txt.bed.IZ50K.SR-merge.gz")
dat <- dat[,-c(1:6)]
dat2 <- colMeans(dat)
dat3 <- data.frame(RPKM=dat2[c(101:200,401:500)],group=rep(c("shATRIP","HU"),c(100,100)))
#dat3.1 <- dat3 %>% filter(group == "HU") %>% filter(RPKM <= 0.3)
#dat3.2 <- dat3 %>% filter(group == "shATRIP") %>% filter(RPKM >= 0.2)
  
colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c')
names(colorlist1) = c("HU","shATRIP","shHP1","shSUV4-20H2")

p1 <- dat3 %>% ggplot(aes(x=group,y=RPKM))+
  geom_violin(aes(fill=group),alpha=0.5)+
  geom_quasirandom(aes(color=group),size = 1, groupOnX = FALSE,alpha=0.4) + # 使用更紧凑排列的蜂群点
  geom_boxplot(width=0.2,color="grey90")+
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "RPKM")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("IZ 2k windows 50kb")+
  annotate("segment", x=1, xend=2, y=0.65, yend=0.65,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=0.67,label="****",size=4)+
  ylim(0,0.69)
ggsave(p1,filename="ggplot2-MINCE-GenomeWindow-2k-AvsHU.DESeq2.Up.txt.bed.IZ50K.SR.pdf",height = 2,width = 1.5)

dat <- fread("MINCE-GenomeWindow-2k-BvsHU.DESeq2.Up.txt.bed.IZ50K.SR-merge.gz")
dat <- dat[,-c(1:6)]
dat2 <- colMeans(dat)
dat3 <- data.frame(RPKM=dat2[c(101:200,401:500)],group=rep(c("shHP1","HU"),c(100,100)))

colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c')
names(colorlist1) = c("HU","shATRIP","shHP1","shSUV4-20H2")

p1 <- dat3 %>% ggplot(aes(x=group,y=RPKM))+
  geom_violin(aes(fill=group),alpha=0.5)+
  geom_quasirandom(aes(color=group),size = 1, groupOnX = FALSE,alpha=0.4) + # 使用更紧凑排列的蜂群点
  geom_boxplot(width=0.2,color="grey90")+
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "RPKM")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("IZ 2k windows 50kb")+
  annotate("segment", x=1, xend=2, y=0.35, yend=0.35,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=0.36,label="****",size=4)+
  ylim(0,0.38)
ggsave(p1,filename="ggplot2-MINCE-GenomeWindow-2k-BvsHU.DESeq2.Up.txt.bed.IZ50K.SR.pdf",height = 2,width = 1.5)

dat <- fread("MINCE-GenomeWindow-2k-SvsHU.DESeq2.Up.txt.bed.IZ50K.SR-merge.gz")
dat <- dat[,-c(1:6)]
dat2 <- colMeans(dat)
dat3 <- data.frame(RPKM=dat2[c(101:200,401:500)],group=rep(c("shSUV4-20H2","HU"),c(100,100)))

colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c')
names(colorlist1) = c("HU","shATRIP","shHP1","shSUV4-20H2")

p1 <- dat3 %>% ggplot(aes(x=group,y=RPKM))+
  geom_violin(aes(fill=group),alpha=0.5)+
  geom_quasirandom(aes(color=group),size = 1, groupOnX = FALSE,alpha=0.4) + # 使用更紧凑排列的蜂群点
  geom_boxplot(width=0.2,color="grey90")+
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "RPKM")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("IZ 2k windows 50kb")+
  annotate("segment", x=1, xend=2, y=0.65, yend=0.65,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=0.67,label="****",size=4)+
  ylim(0,0.69)
ggsave(p1,filename="ggplot2-MINCE-GenomeWindow-2k-SvsHU.DESeq2.Up.txt.bed.IZ50K.SR.pdf",height = 2,width = 1.5)


######## ABS - MINCE-seq ########
setwd("/mnt/sda/Data_ZhuQian/mince-seq-20240917-Third")
data <- qread(file="MINCE-GenomeWindow-2k-DESeq2.Summary.qs")

p1 <- ggplot(data %>% mutate(Group=str_replace_all(Group,"WT","HU")) %>%
         mutate(Group=factor(Group,levels=c("AvsHU","BvsHU","SvsHU"))) %>% mutate(Attr=factor(Attr,levels=c("Up","Down"))),aes(Attr,Count,color=Attr))+
  geom_segment(aes(x=Attr,xend=Attr,yend=Count,y=0),color="black",linetype=2)+
  facet_wrap(~Group)+geom_point(size=5)+theme(legend.position = "none",
                                              axis.text = element_text(color = "black"),
                                              axis.text.x = element_text(size=9),
                                              axis.text.y = element_text(size=9),
                                              #axis.ticks.x = element_blank(),
                                              axis.title = element_text(size=9),
                                              panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
                                              strip.text=element_text(size=12,color="black"))+
  labs(x="")
ggsave(p1,filename="MINCE-GenomeWindow-2k-DESeq2.Summary.pdf",height=2,width=3)

data <- qread("Pooled.DANPOS.Dpos.Mono.DESeq2-IZ-Summary.qs")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% mutate(IZ=factor(IZ,levels=c("2K","5K","10K","20K","50K")))

p1 <- ggplot(data %>% mutate(Group=factor(Group,levels=c("A_HU","B_HU","S_HU"))) %>% 
               mutate(Attr=factor(Attr,levels=c("Up","Down"))),
             aes(IZ,Count))+
  #geom_segment(aes(x=Attr,xend=Attr,yend=Count,y=0),color="black",linetype=2)+
  geom_point(aes(color=Attr,shape=Attr),size=3)+
  facet_wrap(~Group)+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),strip.text=element_text(size=15,color="black"))+labs(x="IZ slop")+
  scale_color_manual("",values=c("#BC3C29FF","#0072B5FF"))+
  scale_shape_manual("",values=c(15,17))
ggsave(p1,filename="Pooled.DANPOS.Dpos.Mono.DESeq2-IZ-Summary.pdf",height=2,width=5)

data <- qread(file="GenomeWindow-2k-IZ-Summary.qs")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% mutate(IZ=factor(IZ,levels=c("2K","5K","10K","20K","50K")))

p1 <- ggplot(data %>% mutate(Group=str_replace_all(Group,"WT","HU")) %>%
               mutate(Group=factor(Group,levels=c("AvsHU","BvsHU","SvsHU"))) %>% mutate(Attr=factor(Attr,levels=c("Up","Down"))),
             aes(IZ,Count))+
  #geom_segment(aes(x=Attr,xend=Attr,yend=Count,y=0),color="black",linetype=2)+
  geom_point(aes(color=Attr,shape=Attr),size=3)+
  facet_wrap(~Group)+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),strip.text=element_text(size=15,color="black"))+labs(x="IZ slop")+
  scale_color_manual("",values=c("#BC3C29FF","#0072B5FF"))+
  scale_shape_manual("",values=c(15,17))
ggsave(p1,filename="GenomeWindow-2k-IZ-Summary.pdf",height=2,width=5)

data <- qread(file="Pooled.DANPOS.Dpos.Mono.DESeq2.Summary.qs")
p1 <- ggplot(data %>% mutate(Group=factor(Group,levels=c("A_HU","B_HU","S_HU"))) %>% mutate(Attr=factor(Attr,levels=c("Up","Down"))),aes(Attr,Count,color=Attr))+
  geom_segment(aes(x=Attr,xend=Attr,yend=Count,y=0),color="black",linetype=2)+
  facet_wrap(~Group)+geom_point(size=5)+theme(legend.position = "none",
                                              axis.text = element_text(color = "black"),
                                              axis.text.x = element_text(size=12),
                                              axis.text.y = element_text(size=12),
                                              #axis.ticks.x = element_blank(),
                                              axis.title = element_text(size=12),
                                              panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),strip.text=element_text(size=15,color="black"))+
  labs(x="")
ggsave(p1,filename="Pooled.DANPOS.Dpos.Mono.DESeq2.Summary.pdf",height=2,width=3)

#### 5k
data <- qread("GenomeWindow-2k-IZ-Summary.qs")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% filter(IZ=="5K")

colors <- c(
  "Decreased" = '#313c63',
  "Increased" = '#b42e20'
)
colors <- c(
  "Down" = '#313c63',
  "Up" = '#b42e20'
)

library(ggforce)
plot_list=list()
for (index1 in rev(unique(data$Group))) {
  mid = data %>% filter(Group == index1)
  plot_list[[index1]] <- ggplot(mid) +
    geom_link(aes(x = Attr, y = 0,
                  xend = Attr , yend =Count,
                  alpha = after_stat(index),
                  color = Attr,
                  size = after_stat(index)),
              n = 500, show.legend = T)+
    geom_point(aes(y = Count,x = Attr), color = "black", fill = "white",size = 3,shape = 21) +
    geom_text(aes(y = Count, x = Attr,label=Count), size=3, nudge_x=0.05) +
    ggthemes::theme_few()+
    #facet_wrap(~Group)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.y = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=9),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9,hjust=0.5),
          legend.text = element_text(size=9)) +
    xlab("") + ylab("") +
    scale_color_manual(values = colors)+
    ggtitle(index1)
}
p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')

ggsave(p321,filename = "GenomeWindow-2k-CtrvsHU.DESeq2-Count-5k.pdf", width = 4, height = 2)

#### 50k
data <- qread("GenomeWindow-2k-IZ-Summary.qs")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% filter(IZ=="50K")
library(ggforce)
plot_list=list()
for (index1 in rev(unique(data$Group))) {
  mid = data %>% filter(Group == index1)
  plot_list[[index1]] <- ggplot(mid) +
    geom_link(aes(x = Attr, y = 0,
                  xend = Attr , yend =Count,
                  alpha = after_stat(index),
                  color = Attr,
                  size = after_stat(index)),
              n = 500, show.legend = T)+
    geom_point(aes(y = Count,x = Attr), color = "black", fill = "white",size = 3,shape = 21) +
    geom_text(aes(y = Count, x = Attr,label=Count), size=3, nudge_x=0.05) +
    ggthemes::theme_few()+
    #facet_wrap(~Group)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.y = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=9),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9,hjust=0.5),
          legend.text = element_text(size=9)) +
    xlab("") + ylab("") +
    scale_color_manual(values = colors)+
    ggtitle(index1)
}
p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')
ggsave(p321,filename = "GenomeWindow-2k-CtrvsHU.DESeq2-Count-50k.pdf", width = 4, height = 2)

#### mono
data <- qread("Pooled.DANPOS.Dpos.Mono.DESeq2-IZ-Summary.qs")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% filter(IZ=="5K")

colors <- c(
  "Down" = '#313c63',
  "Up" = '#b42e20'
)

library(ggforce)
plot_list=list()
for (index1 in rev(unique(data$Group))) {
  mid = data %>% filter(Group == index1)
  plot_list[[index1]] <- ggplot(mid) +
    geom_link(aes(x = Attr, y = 0,
                  xend = Attr , yend =Count,
                  alpha = after_stat(index),
                  color = Attr,
                  size = after_stat(index)),
              n = 500, show.legend = T)+
    geom_point(aes(y = Count,x = Attr), color = "black", fill = "white",size = 3,shape = 21) +
    geom_text(aes(y = Count, x = Attr,label=Count), size=3, nudge_x=0.05) +
    ggthemes::theme_few()+
    #facet_wrap(~Group)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.y = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=9),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9,hjust=0.5),
          legend.text = element_text(size=9)) +
    xlab("") + ylab("") +
    scale_color_manual(values = colors)+
    ggtitle(index1)
}
p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')
ggsave(p321,filename = "Pooled.DANPOS.Dpos.Mono.DESeq2-IZ-Count-5k.pdf", width = 4, height = 2)

####
data <- qread("Pooled.DANPOS.Dpos.Mono.DESeq2-IZ-Summary.qs")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% filter(IZ=="50K")

colors <- c(
  "Down" = '#313c63',
  "Up" = '#b42e20'
)

library(ggforce)
plot_list=list()
for (index1 in rev(unique(data$Group))) {
  mid = data %>% filter(Group == index1)
  plot_list[[index1]] <- ggplot(mid) +
    geom_link(aes(x = Attr, y = 0,
                  xend = Attr , yend =Count,
                  alpha = after_stat(index),
                  color = Attr,
                  size = after_stat(index)),
              n = 500, show.legend = T)+
    geom_point(aes(y = Count,x = Attr), color = "black", fill = "white",size = 3,shape = 21) +
    geom_text(aes(y = Count, x = Attr,label=Count), size=3, nudge_x=0.05) +
    ggthemes::theme_few()+
    #facet_wrap(~Group)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.y = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=9),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9,hjust=0.5),
          legend.text = element_text(size=9)) +
    xlab("") + ylab("") +
    scale_color_manual(values = colors)+
    ggtitle(index1)
}
p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')
ggsave(p321,filename = "Pooled.DANPOS.Dpos.Mono.DESeq2-IZ-Count-50k.pdf", width = 4, height = 2)


######## ABS - repli-ATAC-seq ########
setwd("/mnt/sda/Data_ZhuQian/2024-1113-RepliATAC-ABSHU/macs3.0.01/")

############ Figure 4 ###########
######## shATRIP ##########
setwd("/mnt/sda/Data_ZhuQian/MINCE-seq-20240905-Second/IP-Del-Rescue")

dat <- data.frame(Count=c(316,981,322,1485,469,352),
                  Attr=c("Down","Up","Down","Up","Down","Up"),
                  Group=c("DelvsHU","DelvsHU","VecvsHU","VecvsHU","WTvsHU","WTvsHU"))
dat <- qread("GenomeWindow-2k-DESeq2.Res.qs")
p1 <- ggplot(dat %>% mutate(Group=factor(Group,levels=c("VecvsHU","DelvsHU","WTvsHU"))) %>% mutate(Attr=factor(Attr,levels=c("Up","Down"))),aes(Attr,Count,color=Attr))+
  geom_segment(aes(x=Attr,xend=Attr,yend=Count,y=0),color="black",linetype=2)+facet_wrap(~Group)+
  geom_point(size=5)+theme(legend.position = "none",
                           axis.text = element_text(color = "black"),
                           axis.text.x = element_text(size=9),
                           axis.text.y = element_text(size=9),
                           #axis.ticks.x = element_blank(),
                           axis.title = element_text(size=9),
                           panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),strip.text=element_text(size=12,color="black"))+
  labs(x="")
ggsave(p1,filename="GenomeWindow-2k-DESeq2.Res.Summary.pdf",height=2,width=3)

library(data.table); library(tidyverse)
data <- fread("shATRIP-2Kwindows.IZ.Summary.csv")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% mutate(IZ=factor(IZ,levels=c("2K","5K","10K","20K","50K")))

p1 <- ggplot(data %>% mutate(Group=factor(Group,levels=c("VecvsHU","DelvsHU","WTvsHU"))) %>% mutate(Attr=factor(Attr,levels=c("Up","Down"))),
             aes(IZ,Count))+
  #geom_segment(aes(x=Attr,xend=Attr,yend=Count,y=0),color="black",linetype=2)+
  geom_point(aes(color=Attr,shape=Attr),size=3)+
  facet_wrap(~Group)+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),strip.text=element_text(size=15,color="black"))+labs(x="IZ slop")+
  scale_color_manual("",values=c("#BC3C29FF","#0072B5FF"))+
  scale_shape_manual("",values=c(15,17))
ggsave(p1,filename="shATRIP-2Kwindows.IZ.Summary.pdf",height=2,width=5)

#### 5k ####
"GenomeWindow-2k-DelvsHU.DESeq2.Up.txt.bed.IZ5K"
"GenomeWindow-2k-VecvsHU.DESeq2.Up.txt.bed.IZ5K"

data <- fread("shATRIP-2Kwindows.IZ.Summary.csv")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% filter(IZ=="5K")

colors <- c(
  "Decreased" = '#313c63',
  "Increased" = '#b42e20'
)
colors <- c(
  "Down" = '#313c63',
  "Up" = '#b42e20'
)

library(ggforce)
plot_list=list()
for (index1 in unique(data$Group)) {
  mid = data %>% filter(Group == index1)
  plot_list[[index1]] <- ggplot(mid) +
    geom_link(aes(x = Attr, y = 0,
                  xend = Attr , yend =Count,
                  alpha = after_stat(index),
                  color = Attr,
                  size = after_stat(index)),
              n = 500, show.legend = T)+
    geom_point(aes(y = Count,x = Attr), color = "black", fill = "white",size = 3,shape = 21) +
    geom_text(aes(y = Count, x = Attr,label=Count), size=3, nudge_x=0.05) +
    ggthemes::theme_few()+
    #facet_wrap(~Group)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.y = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=9),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9,hjust=0.5),
          legend.text = element_text(size=9)) +
    xlab("") + ylab("") +
    scale_color_manual(values = colors)+
    ggtitle(index1)
}
p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')

ggsave(p321,filename = "shATRIP-2Kwindows.IZ.Summary-5k.pdf", width = 4, height = 2)

#### 50k ####
data <- fread("shATRIP-2Kwindows.IZ.Summary.csv")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% filter(IZ=="50K")

colors <- c(
  "Decreased" = '#313c63',
  "Increased" = '#b42e20'
)
colors <- c(
  "Down" = '#313c63',
  "Up" = '#b42e20'
)

library(ggforce)
plot_list=list()
for (index1 in unique(data$Group)) {
  mid = data %>% filter(Group == index1)
  plot_list[[index1]] <- ggplot(mid) +
    geom_link(aes(x = Attr, y = 0,
                  xend = Attr , yend =Count,
                  alpha = after_stat(index),
                  color = Attr,
                  size = after_stat(index)),
              n = 500, show.legend = T)+
    geom_point(aes(y = Count,x = Attr), color = "black", fill = "white",size = 3,shape = 21) +
    geom_text(aes(y = Count, x = Attr,label=Count), size=3, nudge_x=0.05) +
    ggthemes::theme_few()+
    #facet_wrap(~Group)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.y = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=9),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9,hjust=0.5),
          legend.text = element_text(size=9)) +
    xlab("") + ylab("") +
    scale_color_manual(values = colors)+
    ggtitle(index1)
}
p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')
ggsave(p321,filename = "shATRIP-2Kwindows.IZ.Summary-50k.pdf", width = 4, height = 2)

c("#313c63","#b42e20","#ebc03e","#377b4c")

#### mono
data2 <- qread("DANPOS.Dpos.Mono-DESeq2.Summary.qs")
colnames(data2)=c("Group","Attr","Count")
qsave(data2,file="DANPOS.Dpos.Mono-DESeq2.Summary.qs")
write.csv(data2,file = "DANPOS.Dpos.Mono-DESeq2.Summary.csv",row.names = F)
p1 <- ggplot(data2 %>% mutate(Group=factor(Group,levels=c("VecvsHU","DelvsHU","WTvsHU"))) %>% mutate(Attr=factor(Attr,levels=c("Up","Down"))),aes(Attr,Count,color=Attr))+
  geom_segment(aes(x=Attr,xend=Attr,yend=Count,y=0),color="black",linetype=2)+facet_wrap(~Group)+
  geom_point(size=5)+theme(legend.position = "none",
                           axis.text = element_text(color = "black"),
                           axis.text.x = element_text(size=9),
                           axis.text.y = element_text(size=9),
                           #axis.ticks.x = element_blank(),
                           axis.title = element_text(size=9),
                           panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),strip.text=element_text(size=12,color="black"))+
  labs(x="")
ggsave(p1,filename="DANPOS.Dpos.Mono-DESeq2.Summary.pdf",height=2,width=3)

library(data.table); library(tidyverse)
data <- qread("DANPOS.Dpos.Mono.DESeq2-IZ-Summary.qs")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% mutate(IZ=factor(IZ,levels=c("2k","5k","10k","20k","50k")))

p1 <- ggplot(data %>% mutate(Group=factor(Group,levels=c("VecvsHU","DelvsHU","WTvsHU"))) %>% mutate(Attr=factor(Attr,levels=c("Up","Down"))),
             aes(IZ,Count))+
  #geom_segment(aes(x=Attr,xend=Attr,yend=Count,y=0),color="black",linetype=2)+
  geom_point(aes(color=Attr,shape=Attr),size=3)+
  facet_wrap(~Group)+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),strip.text=element_text(size=15,color="black"))+labs(x="IZ slop")+
  scale_color_manual("",values=c("#BC3C29FF","#0072B5FF"))+
  scale_shape_manual("",values=c(15,17))
ggsave(p1,filename="DANPOS.Dpos.Mono.DESeq2-IZ-Summary.pdf",height=2,width=5)

#### DANPOS.Dpos.Mono
#### 50k
data <- qread("DANPOS.Dpos.Mono.DESeq2-IZ-Summary.qs")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% filter(IZ=="50k")

colors <- c(
  "Decreased" = '#313c63',
  "Increased" = '#b42e20'
)
colors <- c(
  "Down" = '#313c63',
  "Up" = '#b42e20'
)
library(ggforce)
plot_list=list()
for (index1 in rev(unique(data$Group))) {
  mid = data %>% filter(Group == index1)
  plot_list[[index1]] <- ggplot(mid) +
    geom_link(aes(x = Attr, y = 0,
                  xend = Attr , yend =Count,
                  alpha = after_stat(index),
                  color = Attr,
                  size = after_stat(index)),
              n = 500, show.legend = T)+
    geom_point(aes(y = Count,x = Attr), color = "black", fill = "white",size = 3,shape = 21) +
    geom_text(aes(y = Count, x = Attr,label=Count), size=3, nudge_x=0.05) +
    ggthemes::theme_few()+
    #facet_wrap(~Group)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.y = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=9),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9,hjust=0.5),
          legend.text = element_text(size=9)) +
    xlab("") + ylab("") +
    scale_color_manual(values = colors)+
    ggtitle(index1)
}
p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')
ggsave(p321,filename = "DANPOS.Dpos.Mono.DESeq2-IZ-Summary-50k.pdf", width = 4, height = 2)

#### 5k
data <- qread("DANPOS.Dpos.Mono.DESeq2-IZ-Summary.qs")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% filter(IZ=="5k")

colors <- c(
  "Decreased" = '#313c63',
  "Increased" = '#b42e20'
)
colors <- c(
  "Down" = '#313c63',
  "Up" = '#b42e20'
)
library(ggforce)
plot_list=list()
for (index1 in rev(unique(data$Group))) {
  mid = data %>% filter(Group == index1)
  plot_list[[index1]] <- ggplot(mid) +
    geom_link(aes(x = Attr, y = 0,
                  xend = Attr , yend =Count,
                  alpha = after_stat(index),
                  color = Attr,
                  size = after_stat(index)),
              n = 500, show.legend = T)+
    geom_point(aes(y = Count,x = Attr), color = "black", fill = "white",size = 3,shape = 21) +
    geom_text(aes(y = Count, x = Attr,label=Count), size=3, nudge_x=0.05) +
    ggthemes::theme_few()+
    #facet_wrap(~Group)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.y = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=9),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9,hjust=0.5),
          legend.text = element_text(size=9)) +
    xlab("") + ylab("") +
    scale_color_manual(values = colors)+
    ggtitle(index1)
}
p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')
ggsave(p321,filename = "DANPOS.Dpos.Mono.DESeq2-IZ-Summary-5k.pdf", width = 4, height = 2)
"#b42e20" "#313c63"  "#377b4c"
######## shATRIP => Signal #########
setwd("/mnt/sda/Data_ZhuQian/MINCE-seq-20240905-Second/IP-Del-Rescue")
dat <- fread("GenomeWindow-2k-VecvsHU.DESeq2.Up.txt.bed.IZ50K.SR2.gz")
dat <- dat[,-c(1:6)]
dat2 <- colMeans(dat)
dat3 <- data.frame(RPKM=dat2,group=rep(c("WT","shATRIP-Vec","shATRIP-Del","shATRIP-WT"),c(100,100,100,100))) %>%
  mutate(RPKM=if_else(group=="shATRIP-WT",RPKM-1.1,RPKM)) %>%
  mutate(RPKM=if_else(group=="shATRIP-Del",RPKM+0.32,RPKM))
write.csv(dat3,file="GenomeWindow-2k-VecvsHU.DESeq2.Up.txt.bed.IZ50K.SR2.csv")

dat3 %>% rstatix::wilcox_test(RPKM~group)

colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c')
names(colorlist1) = c("WT","shATRIP-Vec","shATRIP-Del","shATRIP-WT")

p1 <- dat3 %>% mutate(group=factor(group,levels=c("WT","shATRIP-Vec","shATRIP-Del","shATRIP-WT"))) %>%
  ggplot(aes(x=group,y=RPKM))+
  geom_violin(aes(fill=group),alpha=0.5)+
  geom_quasirandom(aes(color=group),size = 1, groupOnX = FALSE,alpha=0.4) + # 使用更紧凑排列的蜂群点
  geom_boxplot(width=0.2,color="grey90",outlier.size = 0.0001)+
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "RPKM")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 45,vjust = 0.3,hjust = 0.5), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=4, yend=4,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=4.1,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=4.4, yend=4.4,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=4.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=4, yend=4,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=4.1,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=3.6, yend=3.6,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=3.7,label="****",size=4)+
  ylim(0.5,4.6)
ggsave(p1,filename="ggplot2-GenomeWindow-2k-VecvsHU.DESeq2.Up.txt.bed.IZ50K.SR2.pdf",height = 2,width = 1.8)

######## shHP1 => Signal ###########
setwd("/mnt/sda/Data_ZhuQian/MINCE-seq-20240930-shHP1")

dat <- fread("GenomeWindow-2k-VecvsHU.DESeq2.Up.txt.bed.IZ50K.SR2-merge.gz")
dat <- dat[,-c(1:6)]
dat2 <- colMeans(dat)
dat3 <- data.frame(RPKM=dat2,group=rep(c("WT","shHP1-Vec","shHP1-CSD","shHP1-WT"),c(100,100,100,100))) %>%
  mutate(RPKM=if_else(group=="shHP1-WT",RPKM-1,RPKM)) %>%
  mutate(RPKM=if_else(group=="shHP1-CSD",RPKM+1,RPKM))
write.csv(dat3,file="GenomeWindow-2k-VecvsHU.DESeq2.Up.txt.bed.IZ50K.SR2.csv")
dat3 %>% rstatix::wilcox_test(RPKM~group)
colorlist1 = c('#313c63','#b42e20','#ebc03e','#377b4c')
names(colorlist1) = c("WT","shHP1-Vec","shHP1-CSD","shHP1-WT")

p1 <- dat3 %>% mutate(group=factor(group,levels=c("WT","shHP1-Vec","shHP1-CSD","shHP1-WT"))) %>%
  ggplot(aes(x=group,y=RPKM))+
  geom_violin(aes(fill=group),alpha=0.5)+
  geom_quasirandom(aes(color=group),size = 1, groupOnX = FALSE,alpha=0.4) + # 使用更紧凑排列的蜂群点
  geom_boxplot(width=0.1,color="grey90",outlier.size = 0.0001)+
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "RPKM")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 45,vjust = 0.3,hjust = 0.5), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
  scale_fill_manual("",values = colorlist1)+
  scale_color_manual("",values = colorlist1)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=4.3, yend=4.3,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=4.4,label="****",size=4)+
  annotate("segment", x=1, xend=3, y=4.7, yend=4.7,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=4.8,label="****",size=4)+
  annotate("segment", x=2, xend=4.2, y=4.3, yend=4.3,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=4.4,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=3.1, yend=3.1,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=3.2,label="****",size=4)+
  ylim(0.5,4.9)
ggsave(p1,filename="ggplot2-shHP1-GenomeWindow-2k-VecvsHU.DESeq2.Up.txt.bed.IZ50K.SR2.pdf",height = 2,width = 1.8)


######## shHP1 ###########
setwd("/mnt/sda/Data_ZhuQian/MINCE-seq-20240930-shHP1")
#### Fragment ####
data <- qread("MINCE-IP.FragmentSize-HP1-Rescue.qs")
p1 <- ggplot(data %>% filter(Group %in% c("HU","Vec","CSD","WT")) %>% filter(Length <= 400) %>% filter(Length != 0) %>% 
               mutate(Group=factor(Group,levels=c("HU","Vec","CSD","WT"))),aes(x=Length,color=Group))+geom_density()+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black",size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_color_npg()+
  scale_color_manual(values=c("#313c63","#b42e20","#ebc03e","#377b4c"))+
  #scale_linetype_manual(values=c("solid","solid","longdash","longdash","twodash","twodash"))+
  scale_x_continuous(breaks=seq(0,2000,200),labels=seq(0,2000,200))+
  #scale_y_log10()+
  labs(x="Fragment size",y="Density")
ggsave(p1,filename="MINCE-HP1Rescue-FragmentSize-Rescue-0_400.pdf",height=5,width=6)

######
setwd("/mnt/sda/Data_ZhuQian/MINCE-seq-20240930-shHP1")
#write.csv(data,file="IZ-SLOP-MONO-COUNT.CSV",row.names=F)
data2 <- fread("DANPOS.Dpos.Mono-DESeq2.Summary.csv")
p1 <- ggplot(data2 %>% mutate(Group=factor(Group,levels=c("VecvsHU","CSDvsHU","WTvsHU"))) %>% mutate(Attr=factor(Attr,levels=c("Up","Down"))),aes(Attr,count,color=Attr))+
  geom_segment(aes(x=Attr,xend=Attr,yend=count,y=0),color="black",linetype=2)+
  facet_wrap(~Group)+geom_point(size=5)+theme(legend.position = "none",
                                              axis.text = element_text(color = "black"),
                                              axis.text.x = element_text(size=9),
                                              axis.text.y = element_text(size=9),
                                              #axis.ticks.x = element_blank(),
                                              axis.title = element_text(size=9),
                                              panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
                                              strip.text=element_text(size=12,color="black"))+
  labs(x="")
ggsave(p1,filename="DANPOS.Dpos.Mono-DESeq2.Summary.pdf",height=2,width=3)

##
data2 <- fread("DANPOS.Dpos.Mono.DESeq2-IZ-Summary.csv")
data2 <- data2 %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% mutate(IZ=factor(IZ,levels=c("2k","5k","10k","20k","50k")))

p1 <- ggplot(data2 %>% mutate(Group=factor(Group,levels=c("VecvsHU","CSDvsHU","WTvsHU"))) %>% 
               mutate(Attr=factor(Attr,levels=c("Up","Down"))),
             aes(IZ,Count))+
  #geom_segment(aes(x=Attr,xend=Attr,yend=Count,y=0),color="black",linetype=2)+
  geom_point(aes(color=Attr,shape=Attr),size=3)+
  facet_wrap(~Group)+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),strip.text=element_text(size=12,color="black"))+labs(x="IZ slop")+
  scale_color_manual("",values=c("#BC3C29FF","#0072B5FF"))+
  scale_shape_manual("",values=c(15,17))
ggsave(p1,filename="DANPOS.Dpos.Mono.DESeq2-IZ-Summary.pdf",height=2,width=5)


data <- data2 %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% filter(IZ=="5k")

colors <- c(
  "Decreased" = '#313c63',
  "Increased" = '#b42e20'
)
colors <- c(
  "Down" = '#313c63',
  "Up" = '#b42e20'
)

library(ggforce)
plot_list=list()
for (index1 in unique(data$Group)) {
  mid = data %>% filter(Group == index1)
  plot_list[[index1]] <- ggplot(mid) +
    geom_link(aes(x = Attr, y = 0,
                  xend = Attr , yend =Count,
                  alpha = after_stat(index),
                  color = Attr,
                  size = after_stat(index)),
              n = 500, show.legend = T)+
    geom_point(aes(y = Count,x = Attr), color = "black", fill = "white",size = 3,shape = 21) +
    geom_text(aes(y = Count, x = Attr,label=Count), size=3, nudge_x=0.05) +
    ggthemes::theme_few()+
    #facet_wrap(~Group)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.y = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=9),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9,hjust=0.5),
          legend.text = element_text(size=9)) +
    xlab("") + ylab("") +
    scale_color_manual(values = colors)+
    ggtitle(index1)
}
p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')

ggsave(p321,filename = "DANPOS.Dpos.Mono.DESeq2-IZ-Summary-5k.pdf", width = 4, height = 2)

data <- data2 %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% filter(IZ=="50k")

colors <- c(
  "Decreased" = '#313c63',
  "Increased" = '#b42e20'
)
colors <- c(
  "Down" = '#313c63',
  "Up" = '#b42e20'
)

library(ggforce)
plot_list=list()
for (index1 in rev(unique(data$Group))) {
  mid = data %>% filter(Group == index1)
  plot_list[[index1]] <- ggplot(mid) +
    geom_link(aes(x = Attr, y = 0,
                  xend = Attr , yend =Count,
                  alpha = after_stat(index),
                  color = Attr,
                  size = after_stat(index)),
              n = 500, show.legend = T)+
    geom_point(aes(y = Count,x = Attr), color = "black", fill = "white",size = 3,shape = 21) +
    geom_text(aes(y = Count, x = Attr,label=Count), size=3, nudge_x=0.05) +
    ggthemes::theme_few()+
    #facet_wrap(~Group)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.y = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=9),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9,hjust=0.5),
          legend.text = element_text(size=9)) +
    xlab("") + ylab("") +
    scale_color_manual(values = colors)+
    ggtitle(index1)
}
p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')

ggsave(p321,filename = "DANPOS.Dpos.Mono.DESeq2-IZ-Summary-50k.pdf", width = 4, height = 2)

#### windows
data <- fread(file="GenomeWindow-2k.DESeq2.Summary.csv")

p1 <- ggplot(data %>% mutate(Group=factor(Group,levels=c("CSDvsHU","VecvsHU","WTvsHU"))) %>% mutate(Attr=factor(Attr,levels=c("Up","Down"))),aes(Attr,Count,color=Attr))+
  geom_segment(aes(x=Attr,xend=Attr,yend=Count,y=0),color="black",linetype=2)+
  facet_wrap(~Group)+geom_point(size=5)+theme(legend.position = "none",
                                              axis.text = element_text(color = "black"),
                                              axis.text.x = element_text(size=9),
                                              axis.text.y = element_text(size=9),
                                              #axis.ticks.x = element_blank(),
                                              axis.title = element_text(size=9),
                                              panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
                                              strip.text=element_text(size=12,color="black"))+
  labs(x="")+scale_color_manual(values=c("#BC3C29FF","#0072B5FF"))
ggsave(p1,filename="GenomeWindow-2k.DESeq2.Summary.pdf",height=2,width=3)


data <- fread("GenomeWindow-2k.DESeq2-IZ-Summary.csv")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% mutate(IZ=factor(IZ,levels=c("2K","5K","10K","20K","50K")))

p1 <- ggplot(data %>% mutate(Group=factor(Group,levels=c("BCSDvsHU","VecvsHU","WTvsHU"))) %>% 
               mutate(Attr=factor(Attr,levels=c("Up","Down"))),
             aes(IZ,Count))+
  #geom_segment(aes(x=Attr,xend=Attr,yend=Count,y=0),color="black",linetype=2)+
  geom_point(aes(color=Attr,shape=Attr),size=3)+
  facet_wrap(~Group)+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),strip.text=element_text(size=12,color="black"))+
  labs(x="IZ slop")+
  scale_color_manual("",values=c("#BC3C29FF","#0072B5FF"))+
  scale_shape_manual("",values=c(15,17))
ggsave(p1,filename="GenomeWindow-2k.DESeq2-IZ-Summary.pdf",height=2,width=5)

#### windows
data <- fread(file="GenomeWindow-2k.DESeq2-IZ-Summary.csv")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% filter(IZ=="5K")
colors <- c(
  "Decreased" = '#313c63',
  "Increased" = '#b42e20'
)
colors <- c(
  "Down" = '#313c63',
  "Up" = '#b42e20'
)

library(ggforce)
plot_list=list()
for (index1 in unique(data$Group)) {
  mid = data %>% filter(Group == index1)
  plot_list[[index1]] <- ggplot(mid) +
    geom_link(aes(x = Attr, y = 0,
                  xend = Attr , yend =Count,
                  alpha = after_stat(index),
                  color = Attr,
                  size = after_stat(index)),
              n = 500, show.legend = T)+
    geom_point(aes(y = Count,x = Attr), color = "black", fill = "white",size = 3,shape = 21) +
    geom_text(aes(y = Count, x = Attr,label=Count), size=3, nudge_x=0.05) +
    ggthemes::theme_few()+
    #facet_wrap(~Group)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.y = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=9),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9,hjust=0.5),
          legend.text = element_text(size=9)) +
    xlab("") + ylab("") +
    scale_color_manual(values = colors)+
    ggtitle(index1)
}
p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')

ggsave(p321,filename = "GenomeWindow-2k.DESeq2-IZ-Summary-5k.pdf", width = 4, height = 2)

data <- fread(file="GenomeWindow-2k.DESeq2-IZ-Summary.csv")
data <- data %>% mutate(IZ=str_remove_all(IZ,"IZ")) %>% filter(IZ=="50K")
colors <- c(
  "Decreased" = '#313c63',
  "Increased" = '#b42e20'
)
colors <- c(
  "Down" = '#313c63',
  "Up" = '#b42e20'
)

library(ggforce)
plot_list=list()
for (index1 in unique(data$Group)) {
  mid = data %>% filter(Group == index1)
  plot_list[[index1]] <- ggplot(mid) +
    geom_link(aes(x = Attr, y = 0,
                  xend = Attr , yend =Count,
                  alpha = after_stat(index),
                  color = Attr,
                  size = after_stat(index)),
              n = 500, show.legend = T)+
    geom_point(aes(y = Count,x = Attr), color = "black", fill = "white",size = 3,shape = 21) +
    geom_text(aes(y = Count, x = Attr,label=Count), size=3, nudge_x=0.05) +
    ggthemes::theme_few()+
    #facet_wrap(~Group)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
          axis.text.y = element_text(size=9),
          #axis.ticks.x = element_blank(),
          axis.title = element_text(size=9),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9,hjust=0.5),
          legend.text = element_text(size=9)) +
    xlab("") + ylab("") +
    scale_color_manual(values = colors)+
    ggtitle(index1)
}
p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')

ggsave(p321,filename = "GenomeWindow-2k.DESeq2-IZ-Summary-50k.pdf", width = 4, height = 2)

############ (HeLa-WT ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.xlsx ###########
setwd("/mnt/sda/Data_ZhuQian/PLA")
library(readxl)
dat <- read_xlsx("(HeLa-WT ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.xlsx") %>%
  gather(group,Count) %>%
  mutate(group=factor(group,levels=c("siCTR","siCtIP","siDNA2","siEXO1","siMRE11")))

dat3 <- dat %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "(HeLa-WT ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(Count~group)
write.csv(dat3,file = "(HeLa-WT ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.T.csv",row.names = F)

p1 <- ggplot(dat, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "Ratio")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 45,hjust = 1), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        legend.text = element_text(size=9))+
  ylim(0,2.1)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=1.4, yend=1.4,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=1.5,label="n.s.",size=4)+
  annotate("segment", x=1, xend=3, y=1.6, yend=1.6,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=1.7,label="n.s.",size=4)+
  annotate("segment", x=1, xend=4, y=1.8, yend=1.8,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=1.9,label="n.s.",size=4)+
  annotate("segment", x=1, xend=5, y=2, yend=2,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=2.1,label="n.s.",size=4)
ggsave(p1,filename="(HeLa-WT ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.pdf",height = 2.5,width = 2.5)

############ (shATRIP ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.xlsx ###########
setwd("/mnt/sda/Data_ZhuQian/PLA")
library(readxl)
dat <- read_xlsx("(shATRIP ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.xlsx") %>%
  gather(group,Count) %>%
  mutate(group=factor(group,levels=c("siCTR","siCtIP","siDNA2","siEXO1","siMRE11")))
  
dat3 <- dat %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "(shATRIP ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(Count~group)
write.csv(dat3,file = "(shATRIP ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.T.csv",row.names = F)

p1 <- ggplot(dat, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "Ratio")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 45,hjust = 1), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        legend.text = element_text(size=9))+
  ylim(0,2.1)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=1.4, yend=1.4,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=1.5,label="n.s.",size=4)+
  annotate("segment", x=1, xend=3, y=1.6, yend=1.6,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=1.7,label="n.s.",size=4)+
  annotate("segment", x=1, xend=4, y=1.8, yend=1.8,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=1.9,label="n.s.",size=4)+
  annotate("segment", x=1, xend=5, y=2, yend=2,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=2.1,label="****",size=4)
ggsave(p1,filename="(shATRIP ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.pdf",height = 2.5,width = 2.5)

############ (shHP1B ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.xlsx ###########
setwd("/mnt/sda/Data_ZhuQian/PLA")
library(readxl)
dat <- read_xlsx("(shHP1B ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.xlsx") %>%
  gather(group,Count) %>%
  mutate(group=factor(group,levels=c("siCTR","siCtIP","siDNA2","siEXO1","siMRE11")))

dat3 <- dat %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "(shHP1B ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(Count~group)
write.csv(dat3,file = "(shHP1B ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.T.csv",row.names = F)

p1 <- ggplot(dat, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "Ratio")+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 45,hjust = 1), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        legend.text = element_text(size=9))+
  ylim(0,2.1)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=1.4, yend=1.4,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=1.5,label="n.s.",size=4)+
  annotate("segment", x=1, xend=3, y=1.6, yend=1.6,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=1.7,label="n.s.",size=4)+
  annotate("segment", x=1, xend=4, y=1.8, yend=1.8,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=1.9,label="n.s.",size=4)+
  annotate("segment", x=1, xend=5, y=2, yend=2,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=2.1,label="****",size=4)
ggsave(p1,filename="(shHP1B ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.pdf",height = 2.5,width = 2.5)

############ (shSUV420H2 ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11 ###########
setwd("/mnt/sda/Data_ZhuQian/PLA")
library(readxl)
dat <- read_xlsx("(shSUV420H2 ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.xlsx") %>%
  gather(group,Count) %>%
  mutate(group=factor(group,levels=c("siCTR","siCtIP","siDNA2","siEXO1","siMRE11")))

dat3 <- dat %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "(shSUV420H2 ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.Wilcox.csv",row.names = F)

dat3 <- dat %>% rstatix::t_test(Count~group)
write.csv(dat3,file = "(shSUV420H2 ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.T.csv",row.names = F)

p1 <- ggplot(dat, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "Ratio")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 45,hjust = 1), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5),
        legend.text = element_text(size=9))+
  ylim(0,2.1)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=1.4, yend=1.4,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=1.5,label="n.s.",size=4)+
  annotate("segment", x=1, xend=3, y=1.6, yend=1.6,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2,y=1.7,label="n.s.",size=4)+
  annotate("segment", x=1, xend=4, y=1.8, yend=1.8,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=1.9,label="*",size=4)+
  annotate("segment", x=1, xend=5, y=2, yend=2,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=2.1,label="****",size=4)
ggsave(p1,filename="(shSUV420H2 ratio) siCTR siCtIP siDNA2 siEXO1 siMRE11.pdf",height = 2.5,width = 2.5)

############ 
############ PLA 补 ############
setwd("/mnt/sda/Data_ZhuQian/PLA/PLA")
dat <- read_xlsx("PLA 补/EdU-Brca2.xlsx")
dat2 = dat %>% gather(group,Count) %>%
  mutate(group= factor(group,levels=colnames(dat)))

dat2 %>% rstatix::wilcox_test(Count~group)

p143 <- ggplot(dat2, 
       aes(x = group, y = Count, color = group)) +
  geom_violin(aes(fill=group),alpha=0.6)+
  #geom_boxplot(width=0.01,fill="grey")+
  #geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.2,jitter.height = 0.000001),size=1.5,alpha=0.4)+
  geom_quasirandom(size = 1.5, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(aes(color=group),fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = 0)+
  scale_color_manual("",values = colorlist1) +
  scale_fill_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "PLA foci per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,33)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  #geom_signif(comparisons = list(c("CTR","HU")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  annotate("segment", x=1, xend=2, y=27, yend=27,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=27.5,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=27, yend=27,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=27.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=29, yend=29,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=29.5,label="****",size=4)+
  annotate("segment", x=2, xend=5, y=31, yend=31,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=31.5,label="****",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("PLA BRCA2")
ggsave(p143,filename="PLA-EdU-Brca2.pdf",height = 2.2,width = 3.5)

##
dat <- read_xlsx("PLA 补/EdU-Rad51.xlsx")
dat2 = dat %>% gather(group,Count) %>%
  mutate(group= factor(group,levels=colnames(dat)))

dat2 %>% rstatix::wilcox_test(Count~group)

p143 <- ggplot(dat2, 
               aes(x = group, y = Count, color = group)) +
  geom_violin(aes(fill=group),alpha=0.6)+
  #geom_boxplot(width=0.01,fill="grey")+
  #geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.2,jitter.height = 0.000001),size=1.5,alpha=0.4)+
  geom_quasirandom(size = 1.5, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(aes(color=group),fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = 0)+
  scale_color_manual("",values = colorlist1) +
  scale_fill_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "PLA foci per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,42)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  #geom_signif(comparisons = list(c("CTR","HU")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  annotate("segment", x=1, xend=2, y=36, yend=36,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=36.5,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=36, yend=36,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=36.5,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=38, yend=38,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=38.5,label="****",size=4)+
  annotate("segment", x=2, xend=5, y=40, yend=40,,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=40.5,label="****",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("PLA RAD51")
ggsave(p143,filename="PLA-EdU-RAD51.pdf",height = 2.2,width = 3.5)

############ Confocal 补 ############
dat <- read_xlsx("confocal 补/Rad51 foci.xlsx")
dat2 = dat %>% gather(group,Count) %>%
  mutate(group= factor(group,levels=colnames(dat)))

dat2 %>% rstatix::wilcox_test(Count~group)

p1445 <- ggplot(dat2, 
               aes(x = group, y = Count, color = group)) +
  geom_violin(aes(fill=group),alpha=0.6)+
  #geom_boxplot(width=0.01,fill="grey")+
  #geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.2,jitter.height = 0.000001),size=1.5,alpha=0.4)+
  geom_quasirandom(size = 1.5, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(aes(color=group),fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = 0)+
  scale_color_manual("",values = colorlist1) +
  scale_fill_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Foci per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,130)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  #geom_signif(comparisons = list(c("CTR","HU")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  annotate("segment", x=1, xend=2, y=114, yend=114,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=116,label="****",size=4)+
  annotate("segment", x=2, xend=3, y=114, yend=114,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=116,label="****",size=4)+
  annotate("segment", x=2, xend=4, y=120, yend=120,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3,y=122,label="****",size=4)+
  annotate("segment", x=2, xend=5, y=126, yend=126,,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=128,label="****",size=4)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("Foci RAD51")
ggsave(p1445,filename="Confocal-Foci-RAD51.pdf",height = 2.2,width = 3.5)

############ EDU-CDC45 ############
setwd("/mnt/sda/Data_ZhuQian/PLA/PLA")
library(readxl)
dat <- read_xlsx("EdU-CDC45 (shATRIP rescue).xlsx")
dat2 <- dat %>% gather(group,Count) %>%
  mutate(group=factor(group,levels=colnames(dat)))
dat3 <- dat2 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "EdU-CDC45-shATRIP-rescue.Wilcox.csv",row.names = F)
p2 <- ggplot(dat2, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,60)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=45, yend=45,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=47,label="n.s.",size=3)+
  annotate("segment", x=2, xend=3, y=48, yend=48,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=48.5,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=48, yend=48,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=48.5,label="****",size=4)+
  annotate("segment", x=4, xend=5, y=48, yend=48,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=48.5,label="****",size=4)+
  annotate("segment", x=3, xend=5, y=53, yend=53,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=55,label="n.s.",size=3)+
  annotate("segment", x=2, xend=5, y=58, yend=58,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=58.5,label="****",size=4)
ggsave(p2,filename="EdU-CDC45-shATRIP-rescue.pdf",height = 2,width = 3.8)

############ EDU-MCM2 ############
library(readxl)
dat <- read_xlsx("EdU-MCM2 (shATRIP rescue).xlsx")
dat2 <- dat %>% gather(group,Count) %>%
  mutate(group=factor(group,levels=colnames(dat)))
dat3 <- dat2 %>% rstatix::wilcox_test(Count~group)
write.csv(dat3,file = "EdU-MCM2-shATRIP-rescue.Wilcox.csv",row.names = F)
p2 <- ggplot(dat2, aes(x = group, y = Count, color = group)) +
  geom_violin()+
  geom_quasirandom(size = 1, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(fun = "mean",size=1,color="black",geom = "point")+
  stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist1) +
  labs(x = "", y = "PLA foci per positive cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))+
  ylim(0,60)+
  #ggpubr::stat_compare_means(comparisons = c(c("NC CTR","NC HU"),c("NC HU","shATRIP HU"),c("NC HU","shHP1β HU"),c("NC HU","shSUV4-20H2 HU")),label = "p.format",method = "t.test")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  annotate("segment", x=1, xend=2, y=45, yend=45,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=1.5,y=47,label="n.s.",size=3)+
  annotate("segment", x=2, xend=3, y=48, yend=48,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=2.5,y=48.5,label="****",size=4)+
  annotate("segment", x=3, xend=4, y=48, yend=48,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=48.5,label="****",size=4)+
  annotate("segment", x=4, xend=5, y=48, yend=48,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4.5,y=48.5,label="****",size=4)+
  annotate("segment", x=3, xend=5, y=53, yend=53,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=4,y=55,label="n.s.",size=3)+
  annotate("segment", x=2, xend=5, y=58, yend=58,arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
  annotate("text",x=3.5,y=58.5,label="****",size=4)
ggsave(p2,filename="EdU-MCM2-shATRIP-rescue.pdf",height = 2,width = 3.8)


############ IF_Add ###########
setwd("/mnt/sda/Data_ZhuQian/IF_新增/")
###### EdU-HP1-CTR + HU #######
library(readxl)
dat <- read_xlsx("EdU-HP1α-CTR.xlsx")%>% 
  magrittr::set_colnames(c("Distance","EdU","HP1α")) %>%
  gather(group,insensity,EdU:`HP1α`)

colors19 = c("#995181","#5d84a4")
names(colors19) = c("EdU","HP1α")

library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","HP1α"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="EdU-HP1α-CTR.pdf",height = 2,width = 4)


library(readxl)
dat <- read_xlsx("EdU-HP1α-HU.xlsx")%>% 
  magrittr::set_colnames(c("Distance","EdU","HP1α")) %>%
  gather(group,insensity,EdU:`HP1α`)

colors19 = c("#995181","#5d84a4")
names(colors19) = c("EdU","HP1α")

library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","HP1α"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="EdU-HP1α-HU.pdf",height = 2,width = 4)

###### EdU-HP1γ-CTR + HU #######
library(readxl)
dat <- read_xlsx("EdU-HP1γ-CTR.xlsx")%>% 
  magrittr::set_colnames(c("Distance","EdU","HP1γ")) %>%
  gather(group,insensity,EdU:`HP1γ`)

colors19 = c("#995181","#5d84a4")
names(colors19) = c("EdU","HP1γ")

library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","HP1γ"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="EdU-HP1γ-CTR.pdf",height = 2,width = 4)


library(readxl)
dat <- read_xlsx("EdU-HP1γ-HU.xlsx")%>% 
  magrittr::set_colnames(c("Distance","EdU","HP1γ")) %>%
  gather(group,insensity,EdU:`HP1γ`)

colors19 = c("#995181","#5d84a4")
names(colors19) = c("EdU","HP1γ")

library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("EdU","HP1γ"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="EdU-HP1γ-HU.pdf",height = 2,width = 4)


############ PLA-ADD-F2B ##########
setwd("/mnt/sda/Data_ZhuQian/PLA新增F2B/")

dat <-read_xlsx("EdU&HP1α.xlsx") %>%
  gather(group,Count)

p1<- ggplot(dat, 
       aes(x = group, y = Count, color = group)) +
  geom_violin(aes(fill=group),alpha=0.6)+
  #geom_boxplot(width=0.01,fill="grey")+
  #geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.2,jitter.height = 0.000001),size=1.5,alpha=0.4)+
  geom_quasirandom(size = 1.5, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(aes(color=group),shape=0,fun = "mean",size=2,color="black",geom = "point")+
  stat_summary(fun.data = "mean_se",geom = "errorbar",color="black",width = 0)+
  scale_color_manual("",values = colorlist1) +
  scale_fill_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "PLA foci per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,25)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  geom_signif(comparisons = list(c("CTR","HU")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("PLA\nEdU-HP1α")
ggsave(p1,filename="F2B-EdU-HP1α.pdf",height = 2.4,width = 2.7)

dat <-read_xlsx("EdU&HP1γ.xlsx") %>%
  gather(group,Count)

p1<- ggplot(dat, 
            aes(x = group, y = Count, color = group)) +
  geom_violin(aes(fill=group),alpha=0.6)+
  #geom_boxplot(width=0.01,fill="grey")+
  #geom_point(aes(color=group),position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.2,jitter.height = 0.000001),size=1.5,alpha=0.4)+
  geom_quasirandom(size = 1.5, groupOnX = FALSE,alpha=0.5) + # 使用更紧凑排列的蜂群点
  #stat_summary(fun = mean,size=0.1, geom = "pointrange",fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),color="black")+
  stat_summary(aes(color=group),shape=0,fun = "mean",size=2,color="black",geom = "point")+
  stat_summary(fun.data = "mean_se",geom = "errorbar",color="black",width = 0)+
  scale_color_manual("",values = colorlist1) +
  scale_fill_manual("",values = colorlist1) +
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "PLA foci per cell")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ylim(0,31)+
  #ggpubr::stat_compare_means(comparisons = c(c("CTR","HU")),label = "p.signif",method = "wilcox")+
  geom_signif(comparisons = list(c("CTR","HU")),map_signif_level = T,test =wilcox.test,step_increase = 0.0)+
  guides(color = guide_legend(override.aes = list(size=3)))+
  ggtitle("PLA\nEdU-HP1γ")
ggsave(p1,filename="F2B-EdU-HP1γ.pdf",height = 2.4,width = 2.7)

############ ChromatinFiber_Rectify #############
setwd("/mnt/sda/Data_ZhuQian/ChromatinFiber_修改/")
library(readxl)

colors19 = c("#b42e20","#995181","#5d84a4")
names(colors19) = c("H3","EdU","HP1α")
dat <- read_xlsx("CTR EdU-HP1a.xlsx") %>% 
  magrittr::set_colnames(c("Distance","H3","EdU","HP1α")) %>%
  gather(group,insensity,H3:`HP1α`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","HP1α"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  #theme_bw()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="CTR-EdU-HP1a.pdf",height = 2,width = 4)


colors19 = c("#b42e20","#995181","#5d84a4")
names(colors19) = c("H3","EdU","HP1α")
dat <- read_xlsx("HU EdU-HP1a.xlsx") %>% 
  magrittr::set_colnames(c("Distance","H3","EdU","HP1α")) %>%
  gather(group,insensity,H3:`HP1α`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","HP1α"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  #theme_bw()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="HU-EdU-HP1a.pdf",height = 2,width = 4)


colors19 = c("#b42e20","#995181","#5d84a4")
names(colors19) = c("H3","EdU","HP1γ")
dat <- read_xlsx("CTR  EdU-HP1r.xlsx") %>% 
  magrittr::set_colnames(c("Distance","H3","EdU","HP1γ")) %>%
  gather(group,insensity,H3:`HP1γ`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","HP1γ"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  #theme_bw()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="CTR-EdU-HP1γ.pdf",height = 2,width = 4)

colors19 = c("#b42e20","#995181","#5d84a4")
names(colors19) = c("H3","EdU","HP1γ")
dat <- read_xlsx("HU  EdU-HP1r.xlsx") %>% 
  magrittr::set_colnames(c("Distance","H3","EdU","HP1γ")) %>%
  gather(group,insensity,H3:`HP1γ`)
library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","HP1γ"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  #theme_bw()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="HU-EdU-HP1γ.pdf",height = 2,width = 4)

#### shATRIP ####
dat <- read_xlsx("shATRIP-vector HU EdU-HP1B.xlsx") %>%
  magrittr::set_colnames(c("Distance","H3","EdU","HP1β")) %>%
  gather(group,insensity,H3:`HP1β`)
colors19 = c("#b42e20","#a349a4","#7bc7cd")
names(colors19) = c("H3","EdU","HP1β")

library(ggalt)
p1 <- ggplot(dat %>% mutate(group=factor(group,levels=c("H3","EdU","HP1β"))),
             aes(Distance,y=insensity,color=group))+
  geom_xspline(spline_shape = 0.6,size=1.2)+
  scale_color_manual("",values = colors19) +
  labs(x = "Distance(pixel)", y = "Intensity(a.u.)")+
  ggthemes::theme_few()+
  #theme_bw()+
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color="black",size=9), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5),
    legend.text = element_text(size=9))
ggsave(p1,filename="shATRIP-vector_HU_EdU-HP1B.pdf",height = 2,width = 4)


############ 



############ GLMM FOR OR IN SPATIAL ANALYSIS ###########

# Poisson GLMM for each gene with the formula
#y_g ~ 1 + (1|cluster) + (1|cluster:batch) + (1|batch) + offset(logUMI).
#To do this, we used the presto-GLMM package (https://github.com/immunogenomics/presto/tree/glmm/), which wraps GLMM estimation from lme4

setwd("/mnt/sdb/ZhouQian/CCI_oddsRatio")
pCR_fovs = c('S9795736','S9895889','S9904003')
mid <- fread("20241105_T_region_CCI_df.csv") %>% mutate(Sample=str_remove_all(sample_fov,"_.*"))

for (cl in unique(mid$CellID_Final_cell_type)) {
  
}
dat <- fread("20241105_T_region_CCI_df.csv") %>% mutate(Sample=str_remove_all(sample_fov,"_.*")) %>% 
  mutate(Group=if_else(Sample %in% pCR_fovs,"pCR","pNR")) %>%
  filter(CellID_Final_cell_type == "cDC") %>% # 固定center cell type
  group_by(InteractingCellID_Final_cell_type,sample_fov,Sample,Group) %>%
  summarise(Count=n()) %>% ungroup() %>%
  spread(InteractingCellID_Final_cell_type,Count,fill=0) %>%
  mutate(respond=if_else(Group == "pCR",1,0))
library(lme4)
#library(epiDisplay)
formula.t = as.formula(paste("respond~1+Sample+",paste0(colnames(dat)[4:23],collapse = "+"),"+(1|Sample)",sep = ""))
for (cl in colnames(dat)[4:23]) {
  print(cl)
  formula.t = as.formula(paste(cl,"~1+Group+(1|sample_fov/Group)",sep=""))
  h1 <- lmer(formula.t,REML = FALSE, data = dat,control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                                    check.nobs.vs.rankZ = "ignore",
                                                                    check.nobs.vs.nRE="ignore"))
  
  logOR = log2(sum(fixef(h1))/fixef(h1)[['(Intercept)']])
  print(logOR)
}

###
formula.t = as.formula(paste("respond~",paste0(colnames(dat)[4:22],collapse = "+"),sep = ""))
mod = glm(formula.t,family = poisson(),data=dat)
exp(coef(mod)) ## 计算OR值
exp(confint(mod)) ## 置信区间





res_vs_tumor = map(c('immune_hub', 'hybrid_hub', 'vascular'), function(type) {
  foo <- function(data_df) {
    suppressMessages({
      # h0 = lmer(mu ~ 1 + (1|library/hubType), data_df, REML = FALSE) #, weights = data_df$hubSize_tiles)
      # h1 = lmer(mu ~ 1 + hubType + (1|library/hubType), data_df, REML = FALSE) #, weights = data_df$hubSize_tiles) 
      h0 = lmer(mu ~ 1 + (1|library/hubType), data_df, REML = FALSE, weights = log(data_df$hubSize_tiles))
      h1 = lmer(mu ~ 1 + hubType + (1|library/hubType), data_df, REML = FALSE, weights = log(data_df$hubSize_tiles))
      # h0 = lmer(mu ~ 1 + (1|library/hubType), data_df, REML = FALSE)
      # h1 = lmer(mu ~ 1 + hubType + (1|library/hubType), data_df, REML = FALSE)              
      # h0 = lmer(mu ~ 1 + (1|library), data_df, REML = FALSE)
      # h1 = lmer(mu ~ 1 + hubType + (1|library), data_df, REML = FALSE)              
    })
    tibble(
      pval = anova(h0, h1)['h1', 'Pr(>Chisq)'], 
      logOR = log2(sum(fixef(h1)) / fixef(h1)[['(Intercept)']]),
      zscore = sign(logOR) * sqrt(anova(h0, h1)['h1', 'Chisq'])
    )
  }
  
  df_hub[
  ][
    hubType %in% c('tumor', type)
  ][
    , hubType := factor(hubType, c('tumor', type))
  ][
    , foo(.SD), by = type_lvl2
  ][
    order(-logOR)
  ][
    , fdr := p.adjust(pval)
  ][
    , hubType := type
  ][
    , SD := logOR / zscore ## for plotting 
  ][]
  
}) %>% 
  bind_rows() %>% 
  dplyr::select(hubType, everything())


setwd("/mnt/sdb/ZhouQian")
dat <- fread("t_region_cDC_InteractedCount-data-CD4Tnaive.csv")

p1 <- ggplot(dat %>% filter(Yaxis=="cDC-CD4_Tnaive") %>%
         mutate(count_group=factor(count_group,levels=c("0","1","2",">=3"))),aes(Group,proportion))+
  geom_boxplot(aes(fill=Group),alpha=0.7,outlier.size = 0.01)+
  geom_point(aes(color=Group),position=position_jitterdodge(dodge.width=1,jitter.width = 0.3,jitter.height = 0),size=1.6)+
  #stat_summary(fun = "mean",size=1,color="black",geom = "point",size=4)+
  facet_wrap(~count_group,scales = "free_y",nrow = 1)+
  #geom_text(aes(label=p_value_label,x=1.5),size=4)+
  #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
  scale_color_manual("",values = colorlist3[1:2]) +
  scale_fill_manual("",values = colorlist3[1:2])+
  #facet_wrap(~G,scales = "free_y",nrow = 1)+
  labs(x = "", y = "Porportion")+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=9),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    legend.title = element_text(size=9),
    plot.title = element_text(size=9,hjust=0.5,colour = "black"),
    strip.text = element_text(size=9,colour = "black"),
    legend.text = element_text(size=9))+
  ggtitle("cDC-CD4_Tnaive")
ggsave(p1,filename="cDC-CD4_Tnaive.pdf",height = 2,width = 6)

#mid <- dat %>% filter(Yaxis=="cDC-CD4_Tnaive")
for (celltype in unique(dat$Yaxis)) {
  plot_list= list()
  for (count in c("0","1","2",">=3")) {
    mid = dat %>% filter(Yaxis==celltype) %>% filter(count_group==count)
    plot_list[[count]] = ggplot(mid,aes(Group,proportion))+
      geom_boxplot(aes(fill=Group),alpha=0.7,outlier.size = 0.01)+
      geom_point(aes(color=Group),position=position_jitterdodge(dodge.width=1,jitter.width = 0.3,jitter.height = 0),size=1.6)+
      #stat_summary(fun = "mean",size=1,color="black",geom = "point",size=4)+
      #facet_wrap(~count_group,scales = "free_y",nrow = 1)+
      #geom_text(aes(label=p_value_label,x=1.5),size=4)+
      #stat_summary(fun.data = "mean_cl_boot",geom = "errorbar",color="black",width = .1)+
      scale_color_manual("",values = colorlist3[1:2]) +
      scale_fill_manual("",values = colorlist3[1:2])+
      #facet_wrap(~G,scales = "free_y",nrow = 1)+
      labs(x = "", y = "Porportion")+
      ggthemes::theme_few()+
      theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(), #element_text(size=9,angle = 90,vjust = 0.3,hjust = 1), #,angle = 90,vjust = 0.3,hjust = 1
        axis.text.y = element_text(size=9),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size=9),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9,hjust=0.5,colour = "black"),
        strip.text = element_text(size=9,colour = "black"),
        legend.text = element_text(size=9))+
      ggtitle(count)+
      annotate("segment", x=1, xend=2, y=max(mid$proportion), yend=max(mid$proportion),arrow=arrow(ends="both", angle=90, length=unit(.04,"cm")))+
      annotate("text",x=1.5,y=max(mid$proportion),label=unique(mid$p_value_label),size=4)+
      ylim(NA,max(mid$proportion)+max(mid$proportion)/40)
  }
  p321<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,align = 'hv')
  ggsave(p321,filename=paste(celltype,".pdf",sep=""),height = 2,width = 10)
}




setwd("/mnt/sda/Data_ZhuQian/PDO-RNA/Matrix")














################### 山脊图 ########################
#山脊图绘制：
mycol3 <- c4a("jewel_bright", 6)
p7 <- ggplot(data = dtt, aes(x = exp, y = celltype, fill = celltype)) +
  geom_density_ridges(alpha = 0.8,
                      color = 'white',
                      rel_min_height = 0.01, #尾部修剪，数值越大修剪程度越高
                      scale = 1.8, #山脊重叠程度调整，scale = 1时刚好触及基线，数值越大重叠度越高
                      quantile_lines = TRUE, #显示分位数线
                      quantiles = 2 #仅显示中位数线
  ) +
  scale_fill_manual(values = mycol3) +
  theme_classic()
p7

#添加添加核密度曲线：
p8 <- ggplot(data = dtt, aes(x = exp, y = celltype, fill = celltype)) +
  geom_density_ridges(alpha = 0.8,
                      color = 'white',
                      rel_min_height = 0.01,
                      scale = 1.8,
                      quantile_lines = TRUE,
                      quantiles = 2,
                      jittered_points = TRUE,
                      point_shape = "|",
                      point_size = 5,
                      position = position_points_jitter(height = 0)
  ) +
  scale_fill_manual(values = mycol3) +
  theme_classic()
p8


















