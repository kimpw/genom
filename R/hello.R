# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

genomic <- function(x,y,z) {
  library(readr)
  library(ggrepel)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(data.table)
  library(magrittr)
  library(tidyverse)
  library(pipeR)
  library(gtools)
  g<-x
  m<-y
  m%<>%rename(gene=ENSG)
  m$CHR<-as.integer(m$CHR)
  g1<-g%>%filter(0<P & P<0.05)
  test2<-g%>%filter(P>0.05)%>%group_by(CHR)%>%sample_frac(0.01)%>%bind_rows(g1)%>%select(-SNP)%>%
    mutate(Tis="Single SNP",
           P=log10(P))%>%filter(CHR!=23)

  tissue<-z
  test<-tissue%>%left_join(m,by='gene')%>%mutate(BP=(START_POS+END_POS)/2,
                                                 P=-log10(pvalue))%>%select(CHR,BP,P,tissue,gene_name)


  df.tmp <- test %>%

    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(test, ., by=c("CHR"="CHR")) %>%

    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot,
            col_cat=ifelse(CHR %%2 ,1,0))%>%
    select(CHR,tissue,P,BPcum,col_cat,gene_name)


  df.tmp2 <- test2 %>%

    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(test2, ., by=c("CHR"="CHR")) %>%

    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot,
            col_cat=ifelse(CHR %%2 ,1,0))%>%
    select(CHR,Tis,P,BPcum,col_cat)

  # Add highlight and annotation information
  #mutate( is_highlight=ifelse(SNP %in% mysnps, "yes", "no")) %>%
  #mutate( is_annotate=ifelse(P < 1e-6, "yes", "no"))

  #df2$CHR<-as.integer(df2$CHR)
  for_tag<-df.tmp%>%filter(P>5)

  axisdf <- df.tmp2 %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


  color_tissue =  c("bisque3","peachpuff2","lightgreen","firebrick4","red","orangered","blue4","darkturquoise","dodgerblue","cyan","lightsteelblue",
                    "lightslateblue","skyblue1", "royalblue","slateblue3","deepskyblue1", "lightpink1","gray32","grey85","goldenrod4","khaki3",
                    "lightgoldenrod","tan2","yellow","tomato","tomato3", "burlywood4","plum","chocolate4","paleturquoise4","palevioletred3","seagreen1",
                    "darkgreen","honeydew4", "cornsilk2","wheat2","sandybrown","lightsalmon3","gold","darkslategrey","forestgreen","lightcoral", "deeppink","gray0")
  # get chromosome center positions for x-axis
  #axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  mypalette <- c("#000000", "#656363")
  labels_cat <- c(sort(unique(as.character(df.tmp$tissue))))

  genom<-ggplot() + theme_bw() +theme(axis.text.x = element_text(size=12, color = 'black'),
                                      axis.text.y = element_text(size=12, color = 'black'),
                                      axis.title.x = element_text(size = 12, face = "bold", color ="black"),
                                      axis.title.y = element_text(size = 12, face = "bold", color ="black"),
                                      axis.ticks.x=element_line(),
                                      legend.key.size = unit(0.5,"cm"),
                                      legend.key = element_rect(size=0.5),
                                      legend.title = element_text(size = 5),
                                      legend.text = element_text(size = 5))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Chromosomes") +
    geom_point(data = df.tmp2, aes(x = BPcum, y = P,size=0.01, color = factor(col_cat))) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    geom_point(data = df.tmp, aes(x = BPcum, y = P,size=0.0005, color = factor(tissue), shape = factor(col_cat)),
               )+
    scale_colour_manual(name = "Tissue Type", values = c("darkgray", "black", color_tissue), labels = c("Single SNP", "Single SNP", labels_cat))+
    guides(shape = "none", size = "none", colour = guide_legend(reverse = TRUE, override.aes = list(size=6))) +

    #geom_label_repel(data = subset(new_data, log_p > gene_tag_p), aes(x = new_pos, y = log_p, label = Gene)) +
    geom_label_repel(data = for_tag, aes(x = BPcum, y = P, label = gene_name)) +
    geom_hline(aes(yintercept = 0), size = 1) +
    geom_hline(yintercept = 5, size = 1) +
    scale_size_area(max_size = 2)
  return(genom)
}
