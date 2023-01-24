
reducereduntGO <- function(df,GO.ID,FDR_col){
  
  require(GO.db)
  ###############################
  ###checking to filter GO:BP where children have better FDR values
  GO_IDs_initial <- data.frame(INITIAL = df[,GO.ID],
                               SUGGESTED = NA)
  ###############################
  message("Filtering step:Filtering results by children terms, evaluating each term.
              Please be patient")
  
  pb <-
    utils::txtProgressBar(min = 0,
                          max = nrow(GO_IDs_initial),
                          style = 3)
  
  SUGGESTED <- list()
  for(m in seq_len(nrow(GO_IDs_initial))){
    #message(m)
    utils::setTxtProgressBar(pb, m)
    #getting children terms per GO:ID enriched
    
    
    children_ids <- stack(as.list(GOBPCHILDREN[GO_IDs_initial$INITIAL[[m]]]))
    #children_ids <-  stack(GOBPCHILDREN_list[GO_IDs_initial$INITIAL[[m]]])
    
    #children terms FDR
    x_fdr_i <- df[df[,GO.ID] %in% children_ids$values,c(GO.ID,FDR_col)]
    colnames(x_fdr_i) <- c("GO.ID","FDR")
    #parental term FDR
    par_fdr <- data.frame(GO.ID = GO_IDs_initial$INITIAL[[m]],
                          FDR = df[,FDR_col][which(df[,GO.ID]==GO_IDs_initial$INITIAL[[m]])])
    
    x_fdr_i <- rbind(x_fdr_i,par_fdr)
    #GO_IDs_initial$SUGGESTED[[m]] <- x_fdr_i$GO.ID[which.min(x_fdr_i$FDR)]
    SUGGESTED[[m]] <- x_fdr_i$GO.ID[which.min(x_fdr_i$FDR)]
    
    rm(children_ids,x_fdr_i,par_fdr)
  };rm(m)
  
  close(pb)
  
  GO_IDs_initial$SUGGESTED <- unlist(SUGGESTED)
  rm(SUGGESTED)
  
  #Getting unique filtered GO terms
  SUGGESTED_GO <- unique(GO_IDs_initial$SUGGESTED)
  #Filtering for all the output and obtain a filter table
  df1 <- df[df[,GO.ID] %in% SUGGESTED_GO,]
  return(df1)
}

# df <- x_s$result
# GO.ID <- "term_id"
# FDR_col <- "p_value"


# x <- reduceredunGO(df,GO.ID,FDR_col)
  