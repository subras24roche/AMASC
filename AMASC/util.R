
match_celltype <- function( cluster_list, cell_type_list ){

    prob_list <- rep(0, length( unique(cluster_list)) )
    for( i in 1:length(unique(cluster_list)) ){
        cluster_this <- which( cluster_list == i)
        prob_list[i]<- length( intersect(cluster_this , cell_type_list ) )/ length(cluster_this)  
    }
    #print(round( prob_list,2) )
    return(prob_list)
} 

label_celltype <- function(protein_data, threshold_protein=7, threshold_cluster=0.35, neighbors=50){

    dm <- t( protein_data )
    Rphenograph_out <- Rphenograph(dm, neighbors)

    feature_names <- rownames(protein_data)
    suffix <- ""
    if( length(grep( "_TotalSeqC", feature_names )) > 0 ) suffix <- "_TotalSeqC"
    if( length(grep( "_TotalSeqB", feature_names )) > 0 ) suffix <- "_TotalSeqB"

    CD3 <- paste0("CD3",suffix)
    CD4 <- paste0("CD4",suffix)  
    CD19 <- paste0("CD19",suffix)  
    CD56 <- paste0("CD56",suffix)  
    CD14 <- paste0("CD14",suffix)  
    CD15 <- paste0("CD15",suffix)  
    CD16 <- paste0("CD16",suffix)  
    CD34 <- paste0("CD34",suffix)
    CD8 <- paste0("CD8a",suffix) 
    if(CD8%in%feature_names==FALSE) CD8 <- paste0("CD8",suffix)

    cd8_list <- Reduce(intersect, list( which( protein_data[CD3, ] > threshold_protein ),
                            which( protein_data[CD8, ] > threshold_protein ),
                            which( protein_data[CD4, ] < threshold_protein )))

    cd4_list <- Reduce(intersect, list( which( protein_data[CD3, ] > threshold_protein ),
                            which( protein_data[CD8, ] < threshold_protein ),
                            which( protein_data[CD4, ] > threshold_protein )))

    dpt_list <- Reduce(intersect, list( which( protein_data[CD3, ] > threshold_protein ),
                            which( protein_data[CD8, ] > threshold_protein ),
                            which( protein_data[CD4, ] > threshold_protein )))

    dnt_list <- Reduce(intersect, list( which( protein_data[CD3, ] > threshold_protein ),
                            which( protein_data[CD8, ] < threshold_protein ),
                            which( protein_data[CD4, ] < threshold_protein )))

    b_list <- which( protein_data[CD19, ] > threshold_protein )

    nk_list <- Reduce(intersect, list( which( protein_data[CD56, ] > threshold_protein ),
                                        which( protein_data[CD3, ] < threshold_protein ) ) )
    if(CD15%in%feature_names==TRUE)
        neu_list <- Reduce(intersect, list( which( protein_data[CD14, ] < threshold_protein ),
                                            which( protein_data[CD15, ] > threshold_protein ),
                                            which( protein_data[CD16, ] > threshold_protein ),
                                            which( protein_data[CD3, ] < threshold_protein ) ) )

    cd14_list <- Reduce(intersect, list( which( protein_data[CD14, ] > threshold_protein ),
                                        which( protein_data[CD3, ] < threshold_protein ) ) )

    dc_list <- Reduce(intersect, list( which( protein_data[CD4, ] > threshold_protein ),
                                    which( protein_data[CD3, ] < threshold_protein ),
                                    which( protein_data[CD14, ] < threshold_protein ) ) )

    if(CD34%in%feature_names==TRUE)
        pro_list <- Reduce(intersect, list( which( protein_data[CD34, ] > threshold_protein ) ) )

    doublet_list <- Reduce(intersect, list( which( protein_data[CD14, ] > threshold_protein ),
                                    which( protein_data[CD3, ] > threshold_protein ) ) )

    ## Annotate the clusters ##

    celltypes <- rep(1, ncol(protein_data))
    names(celltypes) <- colnames(protein_data)
    cluster_list <- membership(Rphenograph_out[[2]])  

    #Other
    celltypes[ which( cluster_list %in% which(match_celltype(cluster_list, doublet_list) > threshold_cluster ) == TRUE ) ] = 1

    #Monocytes
    celltypes[ which( cluster_list %in% which(match_celltype(cluster_list, dc_list) > threshold_cluster ) == TRUE ) ] = 2
    celltypes[ which( cluster_list %in% which(match_celltype(cluster_list, cd14_list) > threshold_cluster ) == TRUE ) ] = 2
    if(CD15%in%feature_names==TRUE)
        celltypes[ which( cluster_list %in% which(match_celltype(cluster_list, neu_list) > threshold_cluster ) == TRUE ) ] = 2

    #B cells
    celltypes[ which( cluster_list %in% which(match_celltype(cluster_list, b_list) > threshold_cluster ) == TRUE ) ] = 4

    #NK cells
    celltypes[ which( cluster_list %in% which(match_celltype(cluster_list, nk_list) > threshold_cluster ) == TRUE ) ] = 5

    #CD8+ T
    celltypes[ which( cluster_list %in% which( match_celltype(cluster_list, cd8_list) > threshold_cluster ) == TRUE ) ] = 6
    celltypes[ which( cluster_list %in% which(match_celltype(cluster_list, dnt_list) > threshold_cluster ) == TRUE ) ] = 6

    #Non-CD8+ T
    celltypes[ which( cluster_list %in% which(match_celltype(cluster_list, cd4_list) > threshold_cluster ) == TRUE ) ] = 8
    celltypes[ which( cluster_list %in% which(match_celltype(cluster_list, dpt_list) > threshold_cluster ) == TRUE ) ] = 8

    #Progenitors
    if(CD34%in%feature_names==TRUE)
        celltypes[ which( cluster_list %in% which(match_celltype(cluster_list, pro_list) > threshold_cluster ) == TRUE ) ] = 16


return(celltypes)
}

#cp<-label_celltype(data.matrix(pbmc_protein), 6,0.4,50)
