
reducereduntGOALL <- function(df,GO.ID,FDR_col,mode){
  
  require(GO.db);require(BBmisc)
  ###############################
  if (is.null(mode) | !mode %in% c("BP","MF","CC")) {
    stop("Please use a valid option")
  }
  ###############################
  if(nrow(df)>0){
    message(paste0("RUNNING MODE GO:",mode))
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
    ###############################
    #obtaining children - parent relationships to test them
    GOBP_stack <- stack(as.list(GOBPCHILDREN[]))
    GOMF_stack <- stack(as.list(GOMFCHILDREN[]))
    GOCC_stack <- stack(as.list(GOCCCHILDREN[]))
                                            
     ###########################################################################
    #if there is not available parent - children relationship,
    #it returns the same GO to avoid errors
    for(m in seq_len(nrow(GO_IDs_initial))){
      #message(m)
      utils::setTxtProgressBar(pb, m)
      #getting children terms per GO:ID enriched fod BP, MF and CC
      
      if(mode=="BP"){
        if(isTRUE(GO_IDs_initial$INITIAL[[m]] %in% GOMF_stack$ind)){
          children_ids <- stack(as.list(GOBPCHILDREN[GO_IDs_initial$INITIAL[[m]]]))
          
        } else{
          children_ids <- data.frame(values=GO_IDs_initial$INITIAL[[m]],
                                     ind=GO_IDs_initial$INITIAL[[m]])
        }
      } else if(mode=="MF"){
        if(isTRUE(GO_IDs_initial$INITIAL[[m]] %in% GOMF_stack$ind)){
          children_ids <- stack(as.list(GOMFCHILDREN[GO_IDs_initial$INITIAL[[m]]]))
          
        } else{
          children_ids <- data.frame(values=GO_IDs_initial$INITIAL[[m]],
                                     ind=GO_IDs_initial$INITIAL[[m]])
        }
      } else if(mode=="CC"){
        if(isTRUE(GO_IDs_initial$INITIAL[[m]] %in% GOCC_stack$ind)){
          children_ids <- stack(as.list(GOCCCHILDREN[GO_IDs_initial$INITIAL[[m]]]))
          
        } else{
          children_ids <- data.frame(values=GO_IDs_initial$INITIAL[[m]],
                                     ind=GO_IDs_initial$INITIAL[[m]])
        }
      }

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
  } else {
    warning("NO RESULTS AVAILABLE")
    df1 <- df
  }
  return(df1)
}

# df <- x_s$result
# GO.ID <- "term_id"
# FDR_col <- "p_value"


# x <- reduceredunGO(df,GO.ID,FDR_col)
