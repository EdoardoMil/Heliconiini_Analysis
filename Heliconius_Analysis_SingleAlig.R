######## ######## ######## ######## ######## ######## ######## 
######## Libraries ######## ######## ######## ######## #######
######## ######## ######## ######## ######## ######## ######## 
library(bio3d)
library(igraph)
library(reshape2)
library(cluster)

######## ######## ######## ######## ######## ######## ######## 
######## FUNCTIONS ######## ######## ######## ######## #######
######## ######## ######## ######## ######## ######## ######## 

Local_Network_Paramters <- function(g){
  
  Betweenness_Centrality <- betweenness(g, v = V(g), directed = FALSE, nobigint = TRUE, normalized = TRUE,  weights = 1/E(g)$weight)
  Closeness_Centrality <- closeness(g, vids = V(g), weights = 1/E(g)$weight, normalized = TRUE)
  Strength_Value <- strength(g)
  Hub_Score_aus <- hub_score(g, weights = E(g)$weight)
  Hub_Score <- Hub_Score_aus$vector
  Clustering_Coefficient <- transitivity(g, type="barrat", weights = E(g)$weight)
  Degree_Value <- degree(g)#, loops = FALSE, normalized = TRUE)
  
  Parameters_DF <- cbind(Betweenness_Centrality, Closeness_Centrality, Strength_Value, Hub_Score, Clustering_Coefficient, Degree_Value)
  Parameters_DF <- as.data.frame(Parameters_DF)
  
  return(Parameters_DF)
  
}

Network_Generation <- function(EdgeList){
  
  EdgeList[,1]=as.character(EdgeList[,1])
  EdgeList[,2]=as.character(EdgeList[,2])
  EdgeList <- as.matrix(EdgeList)
  g <- graph.edgelist(EdgeList[,1:2], directed=FALSE)
  E(g)$weight <- as.numeric(EdgeList[,3])
  
  return(g)
  
}


######## ######## ######## ######## ######## ######## ######## 
######## Importing ######## ######## ######## ######## #######
######## ######## ######## ######## ######## ######## ######## 
dir_input_pdb <- ".."
dir_output_Net <- ".."
dir_output_SeqAlig <- ".."
dir_input_files <- ".."

######## ######## ######## ######## ######## ######## ######## 
########  DATASET  ######## ######## ######## ######## #######
######## ######## ######## ######## ######## ######## ######## 

file_name <- "CocoonaseMainIdTable_mod.csv"
file_dataset <- read.csv(paste(dir_input_files, file_name, sep=""), sep=";")

SubLocus <- unique(as.character(file_dataset$SubLocus))
Locus <- unique(as.character(file_dataset$Locus))
DatasetList <- list()
DatasetList[[1]] <- as.character(file_dataset$Sequence.name)


######## ######## ######## ######## ######## ######## ######## 
########   MAIN    ######## ######## ######## ######## #######
######## ######## ######## ######## ######## ######## ######## 
cutoffDist <- 10
wcutoff <- 0.07
mail_aus <- "edoardo.milanetti@gmail.com"

ControlNames <- c()
l.control <- as.character(file_dataset$Sequence.name)
for(j in 1:length(l.control)){
  pdb_list_aus  <- list.files(paste(dir_input_pdb,l.control[j],"/model",sep=""), pattern = "crderr")
  pdb_aus <- read.pdb(paste(dir_input_pdb,l.control[j],"/model/",pdb_list_aus[1],sep=""))
  if(sum(pdb_aus$calpha) > 190){
    ControlNames <- c(ControlNames,l.control[j])
  } 
}

DatasetList[[1]] <- ControlNames



# Sequence Alignment of PDB Files
for(i in 1:length(DatasetList)){
  round(print(i/length(DatasetList)),3)
  
  listPDB <- list()
  list_aus <- DatasetList[[i]]
  
  df_des_graph_all <- data.frame()
  
  
  if(length(list_aus) >= 2){
    for(j in 1:length(list_aus)){
      pdb_list_aus  <- list.files(paste(dir_input_pdb,list_aus[j],"/model",sep=""), pattern = "crderr")
      pdb_aus <- read.pdb(paste(dir_input_pdb,list_aus[j],"/model/",pdb_list_aus[1],sep=""))
      listPDB[[j]] <- pdb_aus
      
      Net_df <- list()
      for(m in 1:length(pdb_list_aus)){
        pdb_aus <- read.pdb(paste(dir_input_pdb,list_aus[j],"/model/",pdb_list_aus[m],sep=""))
        
        # Distance Matrix
        df_pdb <- pdb_aus$atom
        df_pdb_CA <- df_pdb[df_pdb$elety =="CA", ]
        
          MatDist <- as.matrix(dist(df_pdb_CA[,c("x","y","z")]))
          rownames(MatDist) <- df_pdb_CA$resno 
          colnames(MatDist) <- df_pdb_CA$resno 
          
          # Contact Matrix
          MatDist_Bin <- MatDist
          MatDist_Bin[MatDist_Bin <= cutoffDist] <- 1
          MatDist_Bin[MatDist_Bin > cutoffDist] <- 0
          
          #Graph Theory
          EdgeList_Dist <- melt(MatDist)
          cond1 <- abs(EdgeList_Dist[,1] - EdgeList_Dist[,2]) > 2
          cond2 <- !EdgeList_Dist[,3] == 0
          cond_final <- as.logical(cond1*cond2)
          EdgeList_Dist_final <- EdgeList_Dist[cond_final,]
          EdgeList_Dist_final_inv <-  EdgeList_Dist_final
          EdgeList_Dist_final_inv[,3] <- 1/EdgeList_Dist_final_inv[,3]
          #EdgeList_Dist_final_inv <- EdgeList_Dist_final_inv[EdgeList_Dist_final_inv[,3] > wcutoff,]
          
          
          g_graph <- Network_Generation(EdgeList_Dist_final_inv)
          df_des_graph <- Local_Network_Paramters(g_graph)
          
          Net_df[[m]] <- df_des_graph
        }
        
        Mat1 <- Net_df[[1]]
        for(s in 2:length(Net_df)){
          Mat1 <- Mat1 + Net_df[[s]]
        }
        
        Net_df_ave <- Mat1/length(Net_df)
        df_des_graph$ResNode <- rownames(Net_df_ave)
        df_des_graph$System <- list_aus[j]
        
        df_des_graph_all <- as.matrix(rbind(df_des_graph_all, as.matrix(df_des_graph)))
        
      }
      
      df_des_graph_all <- as.data.frame(df_des_graph_all)
      df_des_graph_all$ResNode <- as.numeric(as.vector(df_des_graph_all$ResNode))
      
      seqAl <- pdbaln(listPDB, web.args = list(email = mail_aus))
      seqAl_res <- seqAl$resno
      seqAl_res[is.na(seqAl_res)] <- 0
      rownames(seqAl_res) <- DatasetList[[i]]
      
      # save file with multiple sequence aligment 
      fileName_seq <- paste("SingleAlig", "_msa", ".csv", sep="")
      write.csv(seqAl_res, paste(dir_output_SeqAlig, fileName_seq), row.names = T)
      
      # Inserting the network parameters into multiple sequence aligment matrix
      descriptors <- colnames(df_des_graph_all)[1:(ncol(df_des_graph_all)-2)]
      
      for(f in 1:length(descriptors)){
        name_des_aus <- descriptors[f]
        cond_col <- as.logical(as.numeric(colnames(df_des_graph_all)==name_des_aus) + as.numeric(colnames(df_des_graph_all)%in%"ResNode") + as.numeric(colnames(df_des_graph_all)%in%"System"))
        df_des_graph_spec <- df_des_graph_all[,cond_col]
        seqAl_des <- matrix(0, nrow = nrow(seqAl_res), ncol = ncol(seqAl_res))
        
        for(r in 1:nrow(seqAl_res)){
          system_aus <- rownames(seqAl_res)[r]
          df_des_graph_spec_aus <- df_des_graph_spec[df_des_graph_spec$System %in% system_aus, ]
          df_des_graph_spec_aus[,1] <- as.numeric(as.vector(df_des_graph_spec_aus[,1]))
          
          for(c in 1:ncol(seqAl_res)){
            res_aus <- seqAl_res[r,c]
            des_aus <- df_des_graph_spec_aus[df_des_graph_spec_aus$ResNode == as.numeric(as.vector(res_aus)),1]
            if(!length(des_aus)==0){
              seqAl_des[r,c] <- des_aus
            }
          }
        }
        
        seqAl_des <- as.data.frame(seqAl_des)
        rownames(seqAl_des) <- rownames(seqAl_res)
        colnames(seqAl_des) <- NULL
        
        # save file with Network parameters
        fileName_des <- paste("SingleAlig_", "_NetDes_", name_des_aus, "_less190.csv", sep="")
        write.csv(seqAl_des, paste(dir_output_Net,fileName_des), row.names = T)
      }
    
  }
}
