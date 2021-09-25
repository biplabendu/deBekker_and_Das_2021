# Enrichment test (GO and PFAM)


# 1. GO enrichment - Function --------------------------------------------------
go_enrichment <- function(geneset,
                          function.dir = ".",
                          data.dir = ".",
                          org = "cflo", 
                          bg = "all", 
                          atleast = 5, 
                          enriched.terms="over") {

  # save the input list of genes for enrichment test
  genes <- geneset

  ## Load the required libraries
  library(tidyverse)

  ## load the selected annotation file
    if (org=="ophio_cflo"){
      
      print("Loading annotation file for Ophiocordyceps camponoti-floridani")
      all_genes <- read.csv(paste0(function.dir,"/functions/func_data/ophio_cflo_annots_robin_ncbi.csv"), 
                            header=T, stringsAsFactors = F, na.strings = c(NA,""," "))
      # define the separator
      separator = "; "
      
      print("Done.")
      
    } else if (org=="ophio_kim"){
      
      print("Loading annotation file for Ophiocordyceps kimflemingae")
      all_genes <- read.csv(paste0(function.dir,"/functions/func_data/ophio_kim_annots_robin_ncbi.csv"), 
                            header=T, stringsAsFactors = F, na.strings = c(NA,""," "))
      # define the separator
      separator = ";"
      
      print("Done.")
      
    } else if (org=="cflo"){
      
      print("Loading annotation file for Camponotus floridanus")
      all_genes <- read.csv(paste0(function.dir,"/functions/func_data/cflo_annots.csv"), 
                            header=T, stringsAsFactors = F, na.strings = c(NA,""," "))
      # define the separator
      separator = "; "
      
      print("Done.")
      
    } else {
      
      print("Invalid option for argument org.")
      print("Select from: ophio_cflo, ophio_kim, cflo")
      stop()
    }
  
  ## make the flattened gene x annotation file
    # Select only the columns that we need
    all_genes <- all_genes[,c("gene_name","GOs")]
    #' Let's replace the NAs in GOs and pfams with "no_annot"
    all_genes[is.na(all_genes)] <- "no_annot"
    #' Let's flatten the files
    all_genes_gos <- all_genes %>%
      dplyr::mutate(go_split = str_split(GOs, separator)) %>%
      unnest() %>%
      #dplyr::select(-GO) %>%
      separate(go_split, c("GO","GO_desc"), sep = "([\\|])", extra = "drop") %>%
      dplyr::select(gene_name, GO, GO_desc) %>%
      unique()
  
  ## define the background geneset to run enrichment against
    if(bg == "all") {
      # save the background data frame with gene_name, GO term, and GO description in an object
      background <- all_genes_gos %>%
        arrange(gene_name)
      
    } else if (bg == "expressed") {
      
      ## If-else statement to load the expressed geneset for the organism
      
      if (org=="ophio_cflo"){
        
        foo <- tbl(dbConnect(RSQLite::SQLite(),paste0(data.dir,"/data/databases/TC6_fungal_data.db")), 
                   "ophio_cflo_expressed_genes") %>% 
                      filter(expressed=="yes") %>% 
                      collect() %>% pull(gene_name) %>% as.character()
        background <- all_genes_gos %>%
          # filter and keep user specified background geneset 
          filter(gene_name %in% foo) %>%
          arrange(gene_name)
      }
      
    } else (
      background <- all_genes_gos %>%
        # filter and keep user specified background geneset 
        filter(gene_name %in% as.character(bg)) %>%
        arrange(gene_name)
    )
  
  ## Need a GO term to GO description file
    go_to_desc <- dplyr::distinct(as.data.frame(all_genes_gos[-1]))

  ## Enrichment to be tested for all GO terms that are:
    ## 1. Present in the test geneset
    ## 2. all unique GO terms each present in at least x number of genes
      annot_terms <- background %>%
        # Keep only the genes in my test geneset
        filter(gene_name %in% genes) %>%
        group_by(GO) %>%
        summarize(num_genes = n()) %>%
        arrange(num_genes) %>%
        # Keep only the GO terms that are annotated in at least 5 genes
        filter(num_genes >= atleast) %>%
        pull(GO)

  # Let's get all the genes in our geneset for each GO term
  # Test geneset dataframe
  df.test <- background %>%
    filter(gene_name %in% genes)
  go_to_genes <- aggregate( .~ GO, df.test, function(x) toString(unique(x)))

  ## Make an empty list that can save your results for each GO term
  df.list <- list()

  ## Print the summary stats for the enrichment test
  print(paste0("Number of genes in background geneset: ", background %>% distinct(gene_name) %>% nrow()))
  print(paste0("Number of genes in the test set: ", length(genes)))
  print("--------------------------------")
  print(paste0("Number of GO terms in background geneset: ", background %>% distinct(GO) %>% nrow()))
  print(paste0("Number of GO terms (at least ", atleast, "genes) in background geneset: ", background %>% group_by(GO) %>% summarise(num_genes = n()) %>% filter(num_genes >= atleast) %>% nrow()))
  print(paste0("Number of GO terms (at least ", atleast, "genes) in test set: ",length(annot_terms)))

  if(length(annot_terms) == 0) {
    return(as.character("There are no GO terms to test enrichment for."))
  } else if (length(annot_terms) >= 1)
  {
    print("Testing for enrichment...")
    # Test the enrichment for each of the GO terms
    for (i in 1:length(annot_terms))
    {
      # get the GO term to be tested for enrichment
      annot <- annot_terms[i]

      # number of DEGs (or genes of interest)
      n_DEG <- length(genes)

      # Number of genes annotated with the GO term in the background gene set
      n_GO <- background %>%
        filter(GO == annot) %>%
        nrow()

      # Number of genes NOT annotated with the GO term in the background gene set
      n_not_GO <- length(unique(background[[1]])) - n_GO

      # Number of genes annotated with the GO term in the test set
      n_GO_DEG <-  df.test %>%
        filter(GO == annot) %>%
        nrow()

      pval <- dhyper(n_GO_DEG,
                     n_GO,
                     n_not_GO,
                     n_DEG, log=F)

      df.list[[i]] <- data.frame(GO = annot_terms[i],
                                 sam_freq = round(n_GO_DEG/n_DEG, 5),
                                 back_freq = round(n_GO/sum(n_GO, n_not_GO), 5),
                                 n_GO_DEG = n_GO_DEG,  # x
                                 n_DEG = n_DEG,      # k
                                 n_GO = n_GO,       # m
                                 #n_not_GO = n_not_GO,  # n
                                 n_all_genes = sum(n_GO,n_not_GO),  # n_all_genes
                                 pVal = pval)


    }

  ## Make the output table:
    df.enriched <- bind_rows(df.list, .id = "column_label") %>%
      dplyr::select(-column_label) %>%
      arrange(pVal) %>%
      mutate(adj_pVal = p.adjust(pVal, "BH")) %>%
      # keeps only the GO terms that are found in the test set
      filter(n_GO_DEG != 0) %>%
      #filter(pVal < 0.1 | adj_pVal < 0.1) %>%   ## I am not filtering anything yet. Return the whole file.
      #left_join(go_to_desc, by="GO") %>%
      left_join(go_to_genes, by="GO") %>%
      mutate(over_under = ifelse(sam_freq > back_freq, "over", "under")) %>%
      dplyr::select(GO, GO_desc, over_under, adj_pVal, everything()) %>%
      arrange(over_under, adj_pVal)

    # if (terms == "over") {
    #   df.enriched <- df.enriched %>%
    #     filter(over_under == "over")
    # }


    return(df.enriched);
  }

  
}

### Usage:
# #Let's test the enichment for some mock geneset
# # Let's load the core dataset to sample some gene names
# load("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_core_datasets.RData")
# mock.genes <- sample(cflo.annots.exp[[1]], 1000)
# mock.enriched.gos <- cflo_go_enrichment(geneset = mock.genes)
# # View the first 5 rows of the results
# mock.enriched.gos[1:5,1:6]

# #Let's try the set of rhythmic genes (RAIN)
# # load datasets and get all the sig. (FDR 5%) rhythmic genes
# load("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_RAIN_datasets.RData")
# for.rhy.genes <- for.rain.atleast4fpkm %>% 
#   mutate(gene_name = rownames(for.rain.atleast4fpkm)) %>% 
#   filter(pVal < 0.05) %>% 
#   pull(gene_name)
# # call the enrichment function
# for.rhy.enriched.gos <- cflo_go_enrichment(geneset = for.rhy.genes)
# for.rhy.enriched.gos %>% 
#   # only keep the over-enriched terms
#   filter(over_under == "over") %>% 
#   # set FDR at 1%
#   filter(adj_pVal < 0.01) %>% 
#   select(1:6)


# 2. plotting GO enrichments ----------------------------------------------

go_enrichment_plot <- function(data, 
                               function.dir = ".",
                               category, 
                               fdr=5, 
                               clean="yes") {

  source(paste0(function.dir,"/functions/theme_publication.R"))

  #Save the data to an object
  df <- data

  # Let's read the file with all Cflo GO terms and their categories
  cflo_gos <- read.csv(paste0(function.dir,"/functions/func_data/distinct_gos_namespace.csv"), header = T, stringsAsFactors = F)
  cflo_gos <- cflo_gos %>%
    dplyr::select(1:3) %>%
    dplyr::select(GO = "GOTerm.identifier", GO_category = "GOTerm.namespace")

  cflo_gos[cflo_gos$GO_category=="biological_process",]$GO_category <- "BP"
  cflo_gos[cflo_gos$GO_category=="cellular_component",]$GO_category <- "CP"
  cflo_gos[cflo_gos$GO_category=="molecular_function",]$GO_category <- "MF"

  col.scheme <- c("#143740",
                  "#286E80",
                  "#3BA4BF",
                  "#4FDBFF",
                  "#47C5E6")


  #Format the dataframe
  df <- df %>%
    # only keep the over-enriched terms
    filter(over_under == "over") %>%
    # remove the NA term
    filter(GO_desc != "NA") %>%
    # set FDR at 1%
    filter(adj_pVal < (fdr/100)) %>%
    # add the rhythmicity scores
    mutate(score = -log(adj_pVal)) %>%
    # add a column containing the GO categories
    left_join(cflo_gos, by="GO") %>%
    # reorder the cols
    dplyr::select(GO, GO_category, everything())


  goplot <- ggplot(df) +
    # set overall appearance of the plot
    theme_Publication() +
    # Define the dependent and independent variables
    aes(x = reorder(GO_desc, score), y = score, fill=GO_category) +
    # Add a line showing the alpha = 0.01 level
    geom_hline(yintercept = -log(0.01), size = 1.2, color = "#C7D7D9", alpha=0.8) +
    # From the defined variables, create a vertical bar chart
    geom_col(position = "dodge", alpha=0.9, size = 1) +
    # Set main and axis titles
    # ggtitle(paste0("GO enrichment | FDR = ",fdr,"%")) +
    xlab("GOs") +
    # add caption
    labs(
      # title = paste0("GO enrichment | FDR = ",fdr,"%"),
      # subtitle = sub,
      caption = "*vertical line denotes FDR of 1%") +
    ylab(expression(-log[10]*' '*q[enriched])) +
    # # add annotations
    # ggrepel::geom_label_repel(aes(label = paste(round((n_GO_DEG/n_GO*100),2), "%", paste(" of ",n_GO,sep=""), sep="")),
    #                           fill = "transparent",
    #                           color = 'black',
    #                           size = 3,
    #                           direction = "x",
    #                           ylim=c(10,max(df$score)),
    #                           point.padding = 0.25,
    #                           label.padding = 0.25,
    #                           segment.color = 'transparent',
    #                           # get rid of the outline for the label
  #                           label.size = NA) +
  theme(legend.position = "bottom") +
    ylim(c(0,max(df$score)+2)) +
    # Add facetting for each GO category
    facet_grid(GO_category ~ ., scales = "free_y", space = "free_y") +
    # Shorten very long labels (GO descriptions)
    scale_x_discrete(label = function(x) stringr::str_trunc(x, 35)) +
    # Flip the x and y axes
    coord_flip() +
    scale_fill_manual(values= setNames(col.scheme, levels(df$GO_category))) +
    # scale_fill_manual(values = setNames(c("lightblue", "darkgreen"), levels(tstat$Hemisphere)))
    theme(strip.background = element_blank(), strip.text = element_blank(), # get rid of facet grid labels
          plot.title = element_text(hjust = 0.5),
          axis.line.y = element_line(colour = "transparent",
                                     size=1),
          legend.title = element_blank(),
          # legend.position = "None",
          plot.caption = element_text(hjust=1),
          axis.title.y = element_blank()) +
    guides(
      fill = guide_legend(
        title = "Legend Title",
        override.aes = aes(label = "")))

  if (clean == "yes") {
    goplot <- goplot
  }

  else(
    goplot <- goplot +
      # add annotations
      ggrepel::geom_label_repel(aes(label = paste(round((n_GO_DEG/n_GO*100),2), "%", paste(" of ",n_GO,sep=""), sep="")),
                                fill = "transparent",
                                color = 'black',
                                size = 3,
                                direction = "x",
                                ylim=c(10,max(df$score)),
                                point.padding = 0.25,
                                label.padding = 0.25,
                                segment.color = 'transparent',
                                # get rid of the outline for the label
                                label.size = NA)
  )


  return(goplot)
}


# # 3. Plotting heat maps ------------------------------------------------------
# 
# cflo_heatmap <- function(geneset, 
#                          cluster.r=F, cluster.c=F,
#                          cutree_cols = 1,
#                          title="Heatmap (z-score)", 
#                          show_rownames=F, show_colnames=F,
#                          annotation =T) {
#   
#   # load libraries
#   library(tidyverse)
#   library(pheatmap)
#   library(viridis)
#   
#   # Let's load the datasets:
#   load(file = "./functions/func_data/TC5_core_datasets.RData")
#   load(file = "./functions/func_data/gene_to_annot.RData")
#   
#   oldnames.for <- names(cflo.zscores.for[,-1])
#   oldnames.nur <- names(cflo.zscores.nur[,-1])
#   newnames.for <- c("2F","4F","6F","8F","10F","12F","14F","16F","18F","20F","22F","24F")
#   newnames.nur <- c("2N","4N","6N","8N","10N","12N","14N","16N","18N","20N","22N","24N")
#   
#   genes = as.character(geneset) %>% unique()
#   
#   ## Load the data - zscores
#   cflo.rhy.exp.for <- cflo.zscores.for %>% 
#     #select(gene_name, X2F:X24F) %>% 
#     #mutate(rownames(cflo.annots.exp) = gene_name) %>% 
#     filter(gene_name %in% genes) %>% 
#     rename_at(vars(oldnames.for), ~ newnames.for)
#   cflo.rhy.exp.nur <- cflo.zscores.nur %>% 
#     #select(gene_name, X2F:X24F) %>% 
#     #mutate(rownames(cflo.annots.exp) = gene_name) %>% 
#     filter(gene_name %in% genes) %>% 
#     rename_at(vars(oldnames.nur), ~ newnames.nur)
#   
#   cflo.exp <- cflo.rhy.exp.for %>% 
#     left_join(cflo.rhy.exp.nur, by="gene_name") %>% 
#     na.omit()
#   
#   # Decide if you want gene name or blast annotation for the gene to be displayed
#   if(annotation == T) {
#     
#     cflo.exp.2 <- cflo.exp %>% left_join(gene_to_annot, by="gene_name")
#     cflo.exp.2 <- 
#       cflo.exp.2 %>% 
#       mutate(annot2 <- paste0(gene_name,"| ",annot))
#     
#     rownames(cflo.exp) = cflo.exp.2$annot2
#     cflo.exp <- data.matrix(cflo.exp[-1])
#     
#   }
#   
#   else if (annotation == F) {
#     rownames(cflo.exp) = cflo.exp$gene_name
#     cflo.exp <- data.matrix(cflo.exp[-1])
#   }
#   
#   
#   # I’ll add some column annotations and create the heatmap.
#   # Annotations for:
#   # 1. Is the sample collected during the light or dark phase? 
#   my_sample_col <- data.frame(caste = rep(c("foragers", "nurses"), c(12,12)),
#                               phase = rep(rep(c("light", "dark", "light"), c(5,6,1)),2))
#   
#   row.names(my_sample_col) <- colnames(cflo.exp)
#   # my_sample_col.nur <- data.frame(phase = rep(c("light", "dark", "light"), c(5,6,1)))
#   # row.names(my_sample_col.nur) <- colnames(cflo.rhy.exp.nur)
#   
#   ## Manual color palette
#   my_colour = list(
#     caste = c(foragers = "#F23030", nurses = "#1A80D9"),
#     phase = c(light = "#F2E205", dark = "#010440")
#     #cluster = c("#D97941", "#F2D5C4", "#BF7C63", "#A6401B"),
#     #overlap = c(no = "white", yes = "#D92818")
#   )
#   # my_colour.nur = list(
#   #   phase = c(light = "#FFF640", dark = "grey60")
#   #   #cluster = c("#D97941", "#F2D5C4", "#BF7C63", "#A6401B"),
#   #   #overlap = c(no = "white", yes = "#D92818")
#   #   )
#   
#   # Color scale
#   my.breaks = seq(min(cflo.exp), max(cflo.exp), by=0.06)
#   
#   heat <- pheatmap(cflo.exp, show_rownames = show_rownames, show_colnames = show_colnames,
#                    #annotation_row = my_gene_col.for[,c("cluster","overlap")], 
#                    annotation_col = my_sample_col,
#                    # cutree_rows = 4,
#                    cutree_cols = cutree_cols,
#                    annotation_colors = my_colour,
#                    border_color=FALSE,
#                    cluster_cols = cluster.c,
#                    # user can specify if clustering occurs or not
#                    cluster_rows = cluster.r,
#                    breaks = my.breaks,
#                    ## color scheme borrowed from: 
#                    color = inferno(length(my.breaks) - 1),
#                    ## remove annotation legend all together
#                    annotation_legend = T,
#                    # treeheight_row = 0, 
#                    # treeheight_col = 0,
#                    # remove the color scale or not
#                    legend = T,
#                    main = as.character(title))
#   
#   # nur.rhy.heat <- pheatmap(cflo.rhy.exp.nur, show_rownames = F, show_colnames = F,
#   #                          # annotation_row = my_gene_col.nur[,c("cluster","overlap")], 
#   #                          annotation_col = my_sample_col.nur,
#   #                          # cutree_rows = 8,
#   #                          # cutree_cols = 2,
#   #                          annotation_colors = my_colour.nur,
#   #                          border_color=FALSE,
#   #                          cluster_cols = F,
#   #                          breaks = my.breaks,
#   #                          ## color scheme borrowed from: 
#   #                          color = inferno(length(my.breaks) - 1),
#   #                          ## remove annotation legend all together/or keep it
#   #                          annotation_legend = T,
#   #                          # treeheight_row = 0, 
#   #                          # treeheight_col = 0,
#   #                          # remove the color scale or not
#   #                          legend = T,
#   #                          main = "Nurses (z-scores)")
#   
#   ## Removing specific elements from the figure
#   # Code borrowed from: https://www.biostars.org/p/351551/
#   
#   # library(grid)
#   # grid.ls(grid.force())
#   # # removing the column annotation 
#   # grid.gedit("col_annotation", gp = gpar(col = "white", fill = "white", text = ""))
#   # # removing the "phase" legend title from the annotation legend
#   # grid.gedit("GRID.text.1225", # THIS NUMBER WILL CHANGE EVERY TIME YOU PLOT
#   #            gp = gpar(col="white"))
#   
#   
#   # arrange the two heatmaps in a single figure
#   # library(gridExtra)
#   # for.nur.rhy.heat <- grid.arrange(grobs = list(for.rhy.heat[[4]], nur.rhy.heat[[4]]), ncol=2)
#   # 
#   # # return the heatmaps as a list
#   # return(list(for.nur.rhy.heat, for.rhy.heat[[4]], nur.rhy.heat[[4]]))
#   
#   return(heat)
# }
# 
# 
# 
