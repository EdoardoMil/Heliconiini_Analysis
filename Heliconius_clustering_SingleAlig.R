######## ######## ######## ######## ######## ######## ######## 
######## Libraries ######## ######## ######## ######## #######
######## ######## ######## ######## ######## ######## ######## 
library(cluster)
library(corrplot)
library(pvclust)
library(rafalib)
library(ggExtra)
library(ggpubr)
library(cowplot)

######## ######## ######## ######## ######## ######## ######## 
######## FUNCTIONS ######## ######## ######## ######## #######
######## ######## ######## ######## ######## ######## ######## 
CountNoZero <- function(x){
  control <- sum(!x==0)
  return(control)
}

MeanReplace <- function(x){
  x_no_zero <- x[x!=0]
  x[x==0] <- mean(x_no_zero)
  return(x)
}

RandomReplace <- function(x){
  x_no_zero <- x[x!=0]
  x[x==0] <- runif(1,min(x_no_zero), max(x_no_zero))
  return(x)
}

######## ######## ######## ######## ######## ######## ######## 
######## Importing ######## ######## ######## ######## #######
######## ######## ######## ######## ######## ######## ######## 
dir_input_net <- "/Users/edo/Documents/Progetti/Heliconius/output_Net_SingleAlig/"
dir_input_files <- "/Users/edo/Documents/Progetti/Heliconius/files/"
dir_output_SeqAlig <- "/Users/edo/Documents/Progetti/Heliconius/output_SeqAlig_SingleAlig/"
NetDes_list <- list.files(dir_input_net)

file_name <- "CocoonaseMainIdTable_mod.csv"
df_dataset <- read.csv(paste(dir_input_files, file_name, sep=""), sep=";")

locusList <- as.character(unique(df_dataset$Locus))
list_Locus_name <- list()
for(i in 1:length(locusList)){
  locusList_aus <- as.character(locusList[i])
  list_Locus_name[[i]] <- df_dataset[as.character(df_dataset$Locus) %in% locusList_aus,"Sequence.name"]
}

file_name_alig <- "SingleAlig_msa.csv"
df_Alig_SingleAlig <- read.csv(paste(dir_output_SeqAlig,file_name_alig, sep=" "), row.names = 1)


######## ######## ######## ######## ######## ######## ######## 
########   MAIN    ######## ######## ######## ######## #######
######## ######## ######## ######## ######## ######## ######## 
df_Alig_SingleAlig_Mat <- as.matrix(df_Alig_SingleAlig)
df_Alig_SingleAlig_Mat <- apply(df_Alig_SingleAlig_Mat,2,as.numeric)


# length of each structure 
df_length <- data.frame()
for(i in 1:nrow(df_Alig_SingleAlig)){
  aus_r <- df_Alig_SingleAlig[i,]
  l_aus <- length(aus_r[!aus_r==0])
  df_length <- as.matrix(rbind(df_length, c(rownames(aus_r),l_aus)))
}



apply(df_Alig_SingleAlig_Mat,1,CountNoZero)
hist(apply(df_Alig_SingleAlig_Mat,1,CountNoZero),100, ylim = c(0,130), xlim=c(50,300), col="gray90", main="Sequence length distribution", xlab="Number of residues", cex.lab=1.4, cex.axis=1.4)
abline(v=190, col="red")

df_Alig_SingleAlig_length <- df_Alig_SingleAlig[!apply(df_Alig_SingleAlig_Mat,1,CountNoZero) < 190,]
df_Alig_SingleAlig_length_Mat <- as.matrix(df_Alig_SingleAlig_length)
df_Alig_SingleAlig_length_Mat <- apply(df_Alig_SingleAlig_length_Mat,2,as.numeric)
hist(apply(df_Alig_SingleAlig_length_Mat,1,CountNoZero),100, ylim = c(0,130), xlim=c(50,300), col="gray90", main="Sequence length distribution", xlab="Number of residues")


NetDes_list <- NetDes_list[grep("190",NetDes_list)]
NetDes_list <- NetDes_list[c(2,3,6)]
NetDes_list <- NetDes_list[c(3,2,1)]

par(mfrow=c(1,3))
ver_name_cluster_small <- c()
for(i in 1:length(NetDes_list)){
  netDes_aus <- NetDes_list[i]
  names_aus <- do.call(rbind,strsplit(NetDes_list,"_"))[,4][i]
  setwd(dir_input_net)
  DesAus <- read.csv(netDes_aus, row.names = 1, header=F)
  
  colnames(DesAus) <- c(paste("V",1:ncol(DesAus),sep=""))
  DesAus_new <- DesAus[rownames(DesAus) %in% df_length[as.numeric(df_length[,2]) > 227 & as.numeric(df_length[,2]) < 245 ,1],]
  DesAus_new <- DesAus_new[,apply(DesAus_new,2,CountNoZero) > 170]
  
  DesAus_new_numeric <- apply(DesAus_new, 2, as.numeric)
  DesAus_new_numeric <- as.data.frame(DesAus_new_numeric)
  rownames(DesAus_new_numeric) <- rownames(DesAus_new)
  DesAus_MeanReplace <- apply(DesAus_new_numeric, 2, MeanReplace)
  
  # corrplot
  mat_corrplot <- DesAus_MeanReplace/max(DesAus_MeanReplace)
  rownames(mat_corrplot) <- NULL
  colnames(mat_corrplot) <- NULL
  #corrplot(mat_corrplot[1:50,1:80], method = "square")
  
  hc_new <- hclust(dist(DesAus_MeanReplace, method = "euclidean"), method = "ward.D2")
  cutree_aus <- cutree(hc_new,2)
  names(cutree_aus[cutree_aus==2])
  ver_name_cluster_small <- c(ver_name_cluster_small, names(cutree_aus[cutree_aus==2]))
  
  condition_color <- rep("gray20", length(rownames(DesAus_MeanReplace)))
  names(condition_color) <- rownames(DesAus_MeanReplace)
  condition_color[names(condition_color) %in% as.character(list_Locus_name[[2]])] <- "red"
  condition_color[names(condition_color) %in% as.character(list_Locus_name[[3]])] <- "blue"
  condition_color[names(condition_color) %in% as.character(list_Locus_name[[4]])] <- "green"
  
  myplclust(hc_new, labels=rownames(DesAus_MeanReplace), lab.col=(as.character(condition_color)), cex.axis=2, cex.lab=2, cex=0.25, hang = 0.005, main = names_aus, mgp=c(2.5,1,0))

  # shilouette
  vet_silh <- c()
  max_cluster <- (nrow(DesAus_MeanReplace)-1)
  dist_aus <- dist(DesAus_MeanReplace, method = "euclidean")
  for(k in 2:50){
    print(k)
    si <- silhouette(cutree(hc_new,k), dist_aus)
    
    summary_si <- summary(si)
    mean_i <- c(summary_si$avg.width)
    vet_silh <- c(vet_silh, mean_i)
    
  }
  
  plot(2:(length(vet_silh)+1), vet_silh, pch=20, cex=3, ylab="Silhouette", xlab="Number of groups", cex.lab=1.5, main = names_aus, ylim=c(0,1))
  lines(2:(length(vet_silh)+1), vet_silh, col="blue", lwd=2)
  
  
  pca_aus <- prcomp((DesAus_MeanReplace))
  #pca_aus <- princomp(DesAus_MeanReplace)
  
  
  Xlim <- c(min(pca_aus$x[,1]) - abs(min(pca_aus$x[,1])*0.15), max(pca_aus$x[,1]) + abs(max(pca_aus$x[,1])*0.15))
  Ylim <- c(min(pca_aus$x[,2]) - abs(min(pca_aus$x[,2])*0.15), max(pca_aus$x[,2]) + abs(max(pca_aus$x[,2])*0.15))

  df_plot <- as.data.frame(cbind(pca_aus$x[,1],pca_aus$x[,2], rep("Coc1A",nrow(pca_aus$x))))
  colnames(df_plot) <- c("PC1","PC2","Sublocus")
  rownames(df_plot) <- rownames(pca_aus$x)
  df_plot$Sublocus <- as.character(df_plot$Sublocus)
  df_plot$PC1 <- as.numeric(as.vector(df_plot$PC1))
  df_plot$PC2 <- as.numeric(as.vector(df_plot$PC2))
  
  df_plot[rownames(df_plot)  %in% as.character(list_Locus_name[[1]]) , "Sublocus"] <- locusList[1]
  df_plot[rownames(df_plot)  %in% as.character(list_Locus_name[[2]]) , "Sublocus"] <- locusList[2]
  df_plot[rownames(df_plot)  %in% as.character(list_Locus_name[[3]]) , "Sublocus"] <- locusList[3]
  df_plot[rownames(df_plot)  %in% as.character(list_Locus_name[[4]]) , "Sublocus"] <- locusList[4]
  
  
  #00316e Coc3
  #00adec Coc2
  #ff4a0d Coc1B
  #f1c49d Coc1A
  
    # Main plot
    pmain <- ggplot(df_plot, aes(x = PC1, y = PC2, color = Sublocus))+
      geom_point(aes(size = 1.5, alpha = 0.3))+
      theme(axis.text=element_text(size=16),
            axis.title=element_text(size=18,face="bold"),
            axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25), 
            panel.background = element_rect(fill = "gray95",
                                            colour = "gray20",
                                            size = 0.5, linetype = "solid"))+
      #ggpubr::color_palette("jco")
      ggpubr::color_palette(c("#f1c49d","#ff4a0d","#00adec","#00316e"))

    pmain
    
    
    # Marginal densities along x axis
    xdens <- axis_canvas(pmain, axis = "x")+
      geom_density(data = df_plot, aes(x = PC1, fill = Sublocus),
                   alpha = 0.7, size = 0.2)+
      #ggpubr::fill_palette("jco")
      ggpubr::fill_palette(c("#f1c49d","#ff4a0d","#00adec","#00316e"))
    
    # Marginal densities along y axis
    # Need to set coord_flip = TRUE, if you plan to use coord_flip()
    ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
      geom_density(data = df_plot, aes(x = PC2, fill = Sublocus),
                   alpha = 0.7, size = 0.2)+
      coord_flip()+
      #ggpubr::fill_palette("jco")
    
    ggpubr::fill_palette(c("#f1c49d","#ff4a0d","#00adec","#00316e"))
    
    p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
    p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
    ggdraw(p2)

  
  # centroids identifications
  
  df_centroids <- data.frame()
  for(k in 1:length(unique(df_plot$Sublocus))){
    sub_aus <- unique(df_plot$Sublocus)[k]
    df_aus <- df_plot[df_plot$Sublocus %in% sub_aus, 1:2]
    centroid_aus <- apply(df_aus,2,mean)
    df_aus_centr <- rbind(df_aus,centroid_aus)
    DistMat <- as.matrix(dist(df_aus_centr))
    name_aus <- names(which.min(DistMat[,nrow(DistMat)][!DistMat[,nrow(DistMat)]==0]))
    vet_info <- c(sub_aus,name_aus)
    df_centroids <- as.matrix(rbind(df_centroids,vet_info))
    
    if(sub_aus == "Coc2"){
      sub_aus_1 <- paste(sub_aus,"_cluster_1", sep="")
      df_aus_1 <- df_aus[df_aus$PC2 > 3,]
      centroid_aus <- apply(df_aus_1,2,mean)
      df_aus_centr <- rbind(df_aus_1,centroid_aus)
      DistMat <- as.matrix(dist(df_aus_centr))
      name_aus <- names(which.min(DistMat[,nrow(DistMat)][!DistMat[,nrow(DistMat)]==0]))
      vet_info <- c(sub_aus_1,name_aus)
      df_centroids <- as.matrix(rbind(df_centroids,vet_info))
      
      sub_aus_2 <- paste(sub_aus,"_cluster_2", sep="")
      df_aus_2 <- df_aus[df_aus$PC2 < 3,]
      centroid_aus <- apply(df_aus_2,2,mean)
      df_aus_centr <- rbind(df_aus_2,centroid_aus)
      DistMat <- as.matrix(dist(df_aus_centr))
      name_aus <- names(which.min(DistMat[,nrow(DistMat)][!DistMat[,nrow(DistMat)]==0]))
      vet_info <- c(sub_aus_2,name_aus)
      df_centroids <- as.matrix(rbind(df_centroids,vet_info))
    }
  }
     
  df_centroids <- as.data.frame(df_centroids)
  rownames(df_centroids) <- NULL
  colnames(df_centroids) <- c("Sublocus","proteinName")
  df_centroids
      
  #loanding...
  load1 <- sort(abs(pca_aus$rotation[,1]), decreasing = T)
  load1[1:10]
  
  load2 <- sort(abs(pca_aus$rotation[,2]), decreasing = T)
  load2[1:10]
  
  load_tot <- sort(abs(pca_aus$rotation[,1])*abs(pca_aus$rotation[,2]), decreasing = T)
  load_tot[1:10]
  
  barplot(load_tot, las=2, cex.axis=1.5, cex.names=0.4, ylim=c(0,0.20))
  load_tot[1:7]
  
  # quanto vale ogni componente...
  ciccio <- summary(pca_aus)
  plot(ciccio$importance[2,], ylab="Proportion of Variance", xlab="Number of components", pch=20, cex=2)
  lines(ciccio$importance[2,], col="gray", lwd=1.8)
 
  plot(ciccio$importance[3,], ylab="Cumulative Proportion", xlab="Number of components", pch=20, cex=2)
  lines(ciccio$importance[3,], col="gray", lwd=1.8)
  
  
  # trova i residui chiave moltiplicanto i loadings della prima e della seconda componente....
  
  
  # Separazione coc2-coc3
  name_coc2 <- as.character(df_dataset[df_dataset$Locus=="Coc2","Sequence.name"])
  name_coc3 <- as.character(df_dataset[df_dataset$Locus=="Coc3","Sequence.name"])
  
  DesAus_MeanReplace_coc2_coc3 <- DesAus_MeanReplace[rownames(DesAus_MeanReplace)  %in% c(name_coc2,name_coc3),]
  pca_coc2_coc3 <- prcomp((DesAus_MeanReplace_coc2_coc3))
  
  df_plot <- as.data.frame(cbind(pca_coc2_coc3$x[,1],pca_coc2_coc3$x[,2], rep("Coc2",nrow(pca_coc2_coc3$x))))
  colnames(df_plot) <- c("PC1","PC2","Sublocus")
  rownames(df_plot) <- rownames(pca_coc2_coc3$x)
  df_plot$Sublocus <- as.character(df_plot$Sublocus)
  df_plot$PC1 <- as.numeric(as.vector(df_plot$PC1))
  df_plot$PC2 <- as.numeric(as.vector(df_plot$PC2))
  
  df_plot[rownames(df_plot)  %in% name_coc2 , "Sublocus"] <- "Coc2"
  df_plot[rownames(df_plot)  %in% name_coc3 , "Sublocus"] <- "Coc3"
  
  library(cowplot) 
  # Main plot
  pmain <- ggplot(df_plot, aes(x = PC1, y = PC2, color = Sublocus))+
    geom_point(aes(size = 1.5, alpha = 0.3))+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18,face="bold"),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25), 
          panel.background = element_rect(fill = "gray95",
                                          colour = "gray20",
                                          size = 0.5, linetype = "solid"))+
  #ggpubr::color_palette("jco")
  ggpubr::color_palette(c("#868686FF","#A73030FF"))
  pmain
  
  # Separazione coc1a-coc1b
  name_coc1a <- as.character(df_dataset[df_dataset$Locus=="Coc1A","Sequence.name"])
  name_coc1b <- as.character(df_dataset[df_dataset$Locus=="Coc1B","Sequence.name"])
  
  DesAus_MeanReplace_coc1a_coc1b <- DesAus_MeanReplace[rownames(DesAus_MeanReplace)  %in% c(name_coc1a,name_coc1b),]
  pca_coc1a_coc1b <- prcomp((DesAus_MeanReplace_coc1a_coc1b))
  
  df_plot <- as.data.frame(cbind(pca_coc1a_coc1b$x[,1],pca_coc1a_coc1b$x[,2], rep("Coc1A",nrow(pca_coc1a_coc1b$x))))
  colnames(df_plot) <- c("PC1","PC2","Sublocus")
  rownames(df_plot) <- rownames(pca_coc1a_coc1b$x)
  df_plot$Sublocus <- as.character(df_plot$Sublocus)
  df_plot$PC1 <- as.numeric(as.vector(df_plot$PC1))
  df_plot$PC2 <- as.numeric(as.vector(df_plot$PC2))
  
  df_plot[rownames(df_plot)  %in% name_coc1a , "Sublocus"] <- "Coc1a"
  df_plot[rownames(df_plot)  %in% name_coc1b , "Sublocus"] <- "Coc1b"
  
  library(cowplot) 
  # Main plot
  pmain <- ggplot(df_plot, aes(x = PC1, y = PC2, color = Sublocus))+
    geom_point(aes(size = 1.5, alpha = 0.3))+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18,face="bold"),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25), 
          panel.background = element_rect(fill = "gray95",
                                          colour = "gray20",
                                          size = 0.5, linetype = "solid"))+
    #ggpubr::color_palette("jco")
    ggpubr::color_palette(c("#0073C2FF","#EFC000FF"))
  pmain
  
  
  }