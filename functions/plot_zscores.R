
# Table of contents -------------------------------------------------------
### Function 1; 
# name: zplot(gene_names, caste); 
# desc: plot zscores 
### Function 2;
# name: multi.plot(plotlist, rows, cols, multipages); 
# desc: plot multiple ggplots into one figure


# Function 1: zplot() -----------------------------------
## name: zplot(gene_names, caste); 
## desc: plot zscores 
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
## 
## gene_names: a character vector; list of names of genes
## caste: a character; possible values include 
##        caste = "for" (plot only forager gene expression),
##        caste = "nur" (plot only nurse gene expression), or
##        caste = "both" (plot both, forager and nurse, gene expression)
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-

zplot <- function(gene_names, caste = "both", lwd=1.5, alpha=0.9) {
  # load the required libraries
  library(tidyverse)
  library(gridExtra)
  library(ggplot2)
  
  # load the core datasets
  load(file = "./functions/func_data/TC5_core_datasets.RData")
  
  # load the theme_Publication() 
  source(file="./functions/theme_publication.R")
  
  #save the current directory 
  current.dir <- getwd()
  #change the directory to the folder where we have our zscores
  #setwd("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data")
  # save the list of gene names in a character vector
  g <- gene_names
  # Read the zscores for caste
  # Read the zscore file
  #foragers_zscores <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data/TC5_zscores_foragers.csv", header = T, stringsAsFactors = F)
  #nurses_zscores <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data/TC5_zscores_nurses.csv", header = T, stringsAsFactors = F)
  foragers_zscores <- cflo.zscores.for
  nurses_zscores <- cflo.zscores.nur
  
  # Read the annotation file for description of the genes
  #all_genes <- read.csv("~/R-scripts_zombie_ant_lab/Enrichment_analysis/Ants/cflo_annotations_expression_v2.csv", header = T, stringsAsFactors = FALSE)
  all_genes <- cflo.annots.exp
  all_genes_annots <- all_genes %>% 
    dplyr::select(gene_name, annot = old_annotation)
  
  # Transforming the data to be able to use ggplot
  # Let's try the logic on a dummy subset
  dummy.for <- foragers_zscores %>%
    filter(gene_name %in% g) %>% 
    gather(ZT, zscore, -1) %>% 
    # arrange(gene_name) %>% 
    arrange(match(gene_name, g)) %>% 
    mutate(caste = "for") %>% 
    mutate(ZT = readr::parse_number(ZT)) %>% 
    left_join(all_genes_annots, by = "gene_name")
  
  dummy.nur <- nurses_zscores %>%
    filter(gene_name %in% g) %>% 
    gather(ZT, zscore, -1) %>% 
    # arrange(gene_name) %>% 
    arrange(match(gene_name, g)) %>% 
    mutate(caste = "nur") %>% 
    mutate(ZT = readr::parse_number(ZT)) %>% 
    left_join(all_genes_annots, by = "gene_name")
  
  # make an if else statement to modify the dataset as per the "caste" parameter
  if (caste == "for" | caste == "forager" | caste == "foragers" | caste == "f"){
    dummy <- dummy.for
    col.scheme <- c("#F23030")
  } 
  else if (caste == "nur" | caste == "nurse" | caste == "nurses" | caste == "n") {
    dummy <- dummy.nur
    col.scheme <- c("#1A80D9")
  } else if (caste == "all" | caste == "both" | caste == "a" | caste == "b") {
    dummy <- rbind(dummy.for, dummy.nur)
    col.scheme <- c("#F23030","#1A80D9")
  } else {
    print("Invalid value for caste. Use one of the following options: for, nur, all.")
    # stop the function and give the user an error.
    stop();
  }
  
  # make the gene_name and caste columns as factors
  dummy[[1]] <- as.factor(dummy[[1]])
  dummy[[4]] <- as.factor(dummy[[4]])
  dummy[[5]] <- as.factor(dummy[[5]])
  
  # Initialize a list to save the plots
  l <- list()
  # Let's plot
  library(ggplot2)
  pd <- position_dodge(0.2)
  l <- lapply(unique(dummy[[1]]), function(i) {
    ggplot(dummy[dummy$gene_name==i,], aes(x=as.numeric(as.character(ZT)), y=zscore)) + 
      #geom_errorbar(aes(ymin=Corrected_Exp-se, ymax=Corrected_Exp+se), width=.1, position=pd) +
      # geom_hline(yintercept=0, color = "red", size = 2) +
      ## if you need highlighting parts of the graph (dark phase in my case)
      geom_rect(aes(xmin = 11.5, xmax = 23.5, ymin = -Inf, ymax = Inf),    
                fill = "lightgrey", alpha = 0.02, color=NA) +
      #facet_grid(~gene_name, scales = "free") +
      facet_wrap(~ gene_name) +
      #ggtitle("TC5 - Z-scores") +
      xlab("") +
      ylab("") +
      theme_Publication() +
      scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
      scale_y_continuous(limits = c(min(dummy[[3]]),max(dummy[[3]]))) +  
      geom_line(data = dummy[dummy$gene_name==i,], position=pd, 
                aes(col=as.factor(caste)), size=lwd, alpha=alpha) +
      # geom_point(position=pd, 
      #            aes(fill=as.factor(caste), shape = as.factor(caste)), size=5, show.legend = F,
      #            color="black", pch=21) +
      #scale_fill_manual(values = c("#F2CB05","#0FBF67")) +
      scale_color_manual(values=col.scheme) + 
      theme(text = element_text(size = 20, colour = 'black'),
            legend.position = "none",
            axis.title.x=element_blank(),
            axis.text.x=element_blank()) +
      # set transparency
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA)) +
      theme(axis.text.y=element_text(color=c("transparent","black","transparent","transparent","black","transparent","transparent","transparent","transparent","black")))
  })
  
  # Let's name the plots with their resp. gene names
  names(l) <- g;
  # let's save the data frame with annotation of the genes into the list as well
  annotations <- all_genes_annots %>% 
    filter(gene_name %in% g) %>% 
    arrange(match(gene_name, g))
  
  # save the list of plots (l) and the annotation table (annotations) into a list
  l2 <- list()
  l2$plots <- l
  l2$annots <- annotations
  
  # set the working directory to the pre-existing one
  setwd(current.dir)
  
  # return the list with all the plots and the data frame containing   
  return(l2);
  
}

# zplot() returns two things:
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
## 
## a list of all the ggplots for each gene ($plots) and 
## a dataframe with the BLAST annotation for each gene ($annots)
## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
## Note: 
## access individual gene plots using: [results]$plots$[gene name]
## access the data frame using: [results]$annots
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-

# Usage for zplot() function
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
###
### get a list of genes that you need to plot; save them in a character vector
#rando.genes <- c("LOC105254558", "LOC105259546", "LOC105256349", "LOC105252017", 
#            "LOC105257043", "LOC105252518", "LOC105250664", "LOC105255790")
### Call the zplot() function
#rando.zplots <- zplot(gene_names = rando.genes)
### Look at the dataframe with the gene annotations
#rando.zplots$annots
### Look at the first gene's expression pattern
#rando.zplots$plots[[1]]
### Look at the same gene's expression pattern by it's name
#rando.zplots$plots$LOC105250664
### For modifying individual ggplots, treat it like an usual ggplot object
# rando.zplots$plots$LOC105250664 +
#   geom_point(aes(col=factor(caste)), size = 4)
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-

# Compiling the multiple plots into one figure
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
## For plotting the multiple gene plots into one or multiple pages, 
## use the multi.plot() function 
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-



# Function 2: multi.plot() -----------------------------------

## name: multi.plot(plotlist, rows, cols, multipages); 
## desc: plot multiple ggplots into one figure 
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
## 
## plotlist = list of ggplots
## rows = number of plots in each column of the page (deafult is 3)
## cols = number of plots in each row of the page (default is 4)
## multipages = binary; 0 indicates print all plots in one page, 1 indicates multi-page plot
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-

multi.plot <- function(plotlist, rows=3, cols=4, multipages="0") {
  # load required packages
  library(gridExtra)
  # coerce the value for multipages into character 
  multi <- as.character(multipages)
  # plot them in 1 page
  if(multipages == "0") {
    # number of total plots?
    tplots <- length(plotlist);
    # number of possible plots given rows and cols?
    pplots <- rows*cols;
    if(tplots <= pplots){
      plots <- do.call(grid.arrange, c(plotlist, ncol=cols, nrow=rows))
    } else {
      print("Can't fit all plots in one page. Change rows or/and cols.")
      stop()
    }
    
  } else if(multipages == "1"){
    plots <- marrangeGrob(plotlist, nrow=rows, ncol=cols)
  } else {
    print("Invalid value for multipages./n 0=one page; 1=multiple pages")
    stop()
  }
  # return the plots object which is of class "arrangelist" "list"
  return(plots)
}

# multi.plot() returns one object of the class "arrangelist" "list"

### Usage for multi.plot() function
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
##
## Call the multi.plot function
## a <- multi.plot(plotlist = rando.plots$plots, rows=2, cols=2, multipages = 1)
## a
## Save the file
## ggsave("name_of_file.png", a, bg = "transparent") # didn't run yet, need to check
## If multipages = 1, you can save the multiple page file into a pdf 
## ggsave("multipage.pdf", a)
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-


exp.plot <- function(gene_names, caste = "both", log=T, lwd=2, alpha=0.9) {
  # load the required libraries
  library(tidyverse)
  library(gridExtra)
  library(ggplot2)
  # load the core datasets
  load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_core_datasets.RData")
  #save the current directory 
  current.dir <- getwd()
  #change the directory to the folder where we have our zscores
  #setwd("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data")
  # save the list of gene names in a character vector
  g <- gene_names
  # Read the zscores for caste
  # Read the zscore file
  #foragers_zscores <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data/TC5_zscores_foragers.csv", header = T, stringsAsFactors = F)
  #nurses_zscores <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data/TC5_zscores_nurses.csv", header = T, stringsAsFactors = F)
  foragers <- cflo.annots.exp %>% dplyr::select(gene_name, X2F:X24F)
  nurses <- cflo.annots.exp %>% dplyr::select(gene_name, X2N:X24N)
  
  # Read the annotation file for description of the genes
  #all_genes <- read.csv("~/R-scripts_zombie_ant_lab/Enrichment_analysis/Ants/cflo_annotations_expression_v2.csv", header = T, stringsAsFactors = FALSE)
  all_genes <- cflo.annots.exp
  all_genes_annots <- all_genes %>% 
    dplyr::select(gene_name, annot = old_annotation)
  
  # Transforming the data to be able to use ggplot
  # Let's try the logic on a dummy subset

  dummy.for <- foragers %>%
    filter(gene_name %in% g) %>% 
    gather(ZT, exp, -1) %>% 
    # mutate(log.exp = log10(exp+1)) %>% 
    # dplyr::select(-exp) %>% 
    arrange(gene_name) %>% 
    mutate(caste = "for") %>% 
    mutate(ZT = readr::parse_number(ZT)) %>% 
    left_join(all_genes_annots, by = "gene_name")
  
  dummy.nur <- nurses %>%
    filter(gene_name %in% g) %>% 
    gather(ZT, exp, -1) %>% 
    # mutate(log.exp = log10(exp+1)) %>% 
    # dplyr::select(-exp) %>% 
    arrange(gene_name) %>% 
    mutate(caste = "nur") %>% 
    mutate(ZT = readr::parse_number(ZT)) %>% 
    left_join(all_genes_annots, by = "gene_name")
  
  # make an if else statement to modify the dataset as per the "caste" parameter
  if (caste == "for" | caste == "forager" | caste == "foragers" | caste == "f"){
    dummy <- dummy.for
    col.scheme <- c("#F23030")
  } 
  else if (caste == "nur" | caste == "nurse" | caste == "nurses" | caste == "n") {
    dummy <- dummy.nur
    col.scheme <- c("#1A80D9")
  } else if (caste == "all" | caste == "both" | caste == "a" | caste == "b") {
    dummy <- rbind(dummy.for, dummy.nur)
    col.scheme <- c("#F23030","#1A80D9")
  } else {
    print("Invalid value for caste. Use one of the following options: for, nur, all.")
    # stop the function and give the user an error.
    stop();
  }
  
  # make the gene_name and caste columns as factors
  dummy[[1]] <- as.factor(dummy[[1]])
  dummy[[4]] <- as.factor(dummy[[4]])
  dummy[[5]] <- as.factor(dummy[[5]])
  
  # Initialize a list to save the plots
  l <- list()
  # Let's plot
  library(ggplot2)
  pd <- position_dodge(0.2)
  
  if(log==T) {
  
  l <- lapply(sort(unique(dummy[[1]])), function(i) {
    ggplot(dummy[dummy$gene_name==i,], aes(x=as.numeric(as.character(ZT)), y=log2(exp+1))) + 
      #geom_errorbar(aes(ymin=Corrected_Exp-se, ymax=Corrected_Exp+se), width=.1, position=pd) +
      # geom_hline(yintercept=0, color = "red", size = 2) +
      ## if you need highlighting parts of the graph (dark phase in my case)
      geom_rect(aes(xmin = 11.5, xmax = 23.5, ymin = -Inf, ymax = Inf),    
                fill = "lightgrey", alpha = 0.02, color=NA) +
      #facet_grid(~gene_name, scales = "free") +
      facet_wrap(~ gene_name, scales = "free_y") +
      #ggtitle("TC5 - Z-scores") +
      xlab("") +
      ylab("") +
      theme_Publication() +
      scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
      # scale_y_continuous(limits = c(min(dummy[[3]]),max(dummy[[3]]))) +  
      geom_line(data = dummy[dummy$gene_name==i,], position=pd, 
                aes(col=as.factor(caste)), size=lwd, alpha=alpha) +
      # geom_point(position=pd, 
      #            aes(fill=as.factor(caste), shape = as.factor(caste)), size=5, show.legend = F,
      #            color="black", pch=21) +
      #scale_fill_manual(values = c("#F2CB05","#0FBF67")) +
      scale_color_manual(values=col.scheme) + 
      theme(text = element_text(size = 15, colour = 'black'),
            legend.position = "none",
            axis.title.x=element_blank(),
            axis.text.x=element_blank()) +
      # set transparency
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA)) +
      theme(axis.text.y=element_text(color=c("transparent","black","transparent","transparent","black","transparent","transparent","transparent","transparent","black")))
  })
  }
  
  else if (log==F) {
    l <- lapply(sort(unique(dummy[[1]])), function(i) {
      ggplot(dummy[dummy$gene_name==i,], aes(x=as.numeric(as.character(ZT)), y=exp)) + 
        #geom_errorbar(aes(ymin=Corrected_Exp-se, ymax=Corrected_Exp+se), width=.1, position=pd) +
        # geom_hline(yintercept=0, color = "red", size = 2) +
        ## if you need highlighting parts of the graph (dark phase in my case)
        geom_rect(aes(xmin = 11.5, xmax = 23.5, ymin = -Inf, ymax = Inf),    
                  fill = "lightgrey", alpha = 0.02, color=NA) +
        #facet_grid(~gene_name, scales = "free") +
        facet_wrap(~ gene_name, scales = "free_y") +
        #ggtitle("TC5 - Z-scores") +
        xlab("") +
        ylab("") +
        theme_Publication() +
        scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
        # scale_y_continuous(limits = c(min(dummy[[3]]),max(dummy[[3]]))) +  
        geom_line(data = dummy[dummy$gene_name==i,], position=pd, 
                  aes(col=as.factor(caste)), size=lwd, alpha=alpha) +
        # geom_point(position=pd, 
        #            aes(fill=as.factor(caste), shape = as.factor(caste)), size=5, show.legend = F,
        #            color="black", pch=21) +
        #scale_fill_manual(values = c("#F2CB05","#0FBF67")) +
        scale_color_manual(values=col.scheme) + 
        theme(text = element_text(size = 15, colour = 'black'),
              legend.position = "none",
              axis.title.x=element_blank(),
              axis.text.x=element_blank()) +
        # set transparency
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA)) +
        theme(axis.text.y=element_text(color=c("black","transparent","black","transparent","transparent","transparent","transparent","black")))
    })
  }
  
  # Let's name the plots with their resp. gene names
  names(l) <- sort(g);
  # let's save the data frame with annotation of the genes into the list as well
  annotations <- all_genes_annots %>% 
    filter(gene_name %in% g)
  
  # save the list of plots (l) and the annotation table (annotations) into a list
  l2 <- list()
  l2$plots <- l
  l2$annots <- annotations
  
  # set the working directory to the pre-existing one
  setwd(current.dir)
  
  # return the list with all the plots and the data frame containing   
  return(l2);
  
}

# zplot() returns two things:
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
## 
## a list of all the ggplots for each gene ($plots) and 
## a dataframe with the BLAST annotation for each gene ($annots)
## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
## Note: 
## access individual gene plots using: [results]$plots$[gene name]
## access the data frame using: [results]$annots
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-

# Usage for zplot() function
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
###
### get a list of genes that you need to plot; save them in a character vector
#rando.genes <- c("LOC105254558", "LOC105259546", "LOC105256349", "LOC105252017", 
#            "LOC105257043", "LOC105252518", "LOC105250664", "LOC105255790")
### Call the zplot() function
#rando.zplots <- zplot(gene_names = rando.genes)
### Look at the dataframe with the gene annotations
#rando.zplots$annots
### Look at the first gene's expression pattern
#rando.zplots$plots[[1]]
### Look at the same gene's expression pattern by it's name
#rando.zplots$plots$LOC105250664
### For modifying individual ggplots, treat it like an usual ggplot object
# rando.zplots$plots$LOC105250664 +
#   geom_point(aes(col=factor(caste)), size = 4)
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-

# Compiling the multiple plots into one figure
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
## For plotting the multiple gene plots into one or multiple pages, 
## use the multi.plot() function 
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-


# Function name: stacked.zplot
# USage: to plot multiple diel gene expression (of Z scores) stacked on top of each other;
# Application: could be used to report similar diel expression patterns for multiple genes

stacked.zplot <- 
  function(gene_names, caste = "both", lwd=1.5, alpha=0.75, bg.alpha = 0.1) {
    
    # load the required libraries
    library(tidyverse)
    library(gridExtra)
    library(ggplot2)
    
    # load the core datasets
    load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_core_datasets.RData")
    
    #save the current directory 
    current.dir <- getwd()
    #change the directory to the folder where we have our zscores
    #setwd("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data")
    # save the list of gene names in a character vector
    g <- gene_names
    
    # if (log=F) {
    # Read the zscores for caste
    # Read the zscore file
    #foragers_zscores <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data/TC5_zscores_foragers.csv", header = T, stringsAsFactors = F)
    #nurses_zscores <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data/TC5_zscores_nurses.csv", header = T, stringsAsFactors = F)
    foragers_zscores <- cflo.zscores.for
    nurses_zscores <- cflo.zscores.nur
    # }
    
    # else {
    #   foragers_zscores <-  cflo.annots.exp %>% dplyr::select(gene_name, X2F:X24F)
    #   foragers_zscores[-1] <- log2(foragers_zscores[-1] + 1)
    #   nurses_zscores <-  cflo.annots.exp %>% dplyr::select(gene_name, X2N:X24N)
    #   nurses_zscores[-1] <- log2(nurses_zscores[-1] + 1)
    # }
    
    # Read the annotation file for description of the genes
    #all_genes <- read.csv("~/R-scripts_zombie_ant_lab/Enrichment_analysis/Ants/cflo_annotations_expression_v2.csv", header = T, stringsAsFactors = FALSE)
    all_genes <- cflo.annots.exp
    all_genes_annots <- all_genes %>% 
      dplyr::select(gene_name, annot = old_annotation)
    
    # Transforming the data to be able to use ggplot
    # Let's try the logic on a dummy subset
    dummy.for <- foragers_zscores %>%
      filter(gene_name %in% g) %>% 
      gather(ZT, zscore, -1) %>% 
      arrange(match(gene_name, g)) %>% 
      mutate(caste = "for") %>% 
      mutate(ZT = readr::parse_number(ZT)) %>% 
      left_join(all_genes_annots, by = "gene_name")
    
    dummy.nur <- nurses_zscores %>%
      filter(gene_name %in% g) %>% 
      gather(ZT, zscore, -1) %>% 
      arrange(match(gene_name, g)) %>% 
      mutate(caste = "nur") %>% 
      mutate(ZT = readr::parse_number(ZT)) %>% 
      left_join(all_genes_annots, by = "gene_name")
    
    # make an if else statement to modify the dataset as per the "caste" parameter
    if (caste == "for" | caste == "forager" | caste == "foragers" | caste == "f"){
      dummy <- dummy.for
      col.scheme <- c("#D91A2A")
    } 
    else if (caste == "nur" | caste == "nurse" | caste == "nurses" | caste == "n") {
      dummy <- dummy.nur
      col.scheme <- c("#2182BF")
    } else if (caste == "all" | caste == "both" | caste == "a" | caste == "b") {
      dummy <- rbind(dummy.for, dummy.nur)
      col.scheme <- c("#D91A2A","#2182BF")
    } else {
      print("Invalid value for caste. Use one of the following options: for, nur, all.")
      # stop the function and give the user an error.
      stop();
    }
    
    # make the gene_name and caste columns as factors
    dummy[[1]] <- as.factor(dummy[[1]]) # gene name
    dummy[[4]] <- as.factor(dummy[[4]]) # caste
    dummy[[5]] <- as.factor(dummy[[5]]) # annotation column
    
    dummy.summary <- 
      dummy %>% 
      group_by(caste, ZT) %>% 
      summarise(mean_zscore = mean(zscore))
    
    # Initialize a list to save the plots
    l <- list()
    # Let's plot
    library(ggplot2)
    pd <- position_dodge(0.2)
    l <- lapply(sort(unique(dummy[[4]])), function(i) {
      
      # Define the dataset
      ggplot() + 
        
        # indicate light-dark phase
        geom_rect(aes(xmin = 11.5, xmax = 23.5, ymin = -Inf, ymax = Inf),
                  fill = "lightgrey", alpha = 0.3, color=NA) +
        
        # Plot the individual gene expressions (zscores)
        geom_line(data=dummy[dummy$caste==i,], 
                  aes(x=as.numeric(as.character(ZT)), y=zscore, group = gene_name),
                  size=0.5, alpha=bg.alpha, 
                  col = ifelse(i=="for",
                                  c(col.scheme[[1]], col.scheme[[2]]),
                                  c(col.scheme[[2]], col.scheme[[1]]))) +
        
        # Plot the mean gene expression for the given gene list
        geom_line(data = dummy.summary[dummy.summary$caste == i,],
                  aes(x=as.numeric(as.character(ZT)), y=mean_zscore), 
                  col = "#3C3C40",
                  size = lwd, alpha = alpha) +
        
        
        # Set the theme
        xlab("") +
        ylab("") +
        theme_Publication() +
        scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
        scale_y_continuous(limits = c(min(dummy[[3]]),max(dummy[[3]]))) + 
        theme(text = element_text(size = 20, colour = 'black'),
              legend.position = "none",
              axis.title.x=element_blank(),
              axis.text.x=element_blank()) +
          # set transparency
        theme( 
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA)) 
        # show only a subset of the y-axis tick labels
        # theme(axis.text.y=element_text(color=c("transparent","black","transparent",
        #                                        "transparent","black","transparent",
        #                                        "transparent","transparent","transparent",
        #                                        "black")))
    })
    
    # Let's name the plots with their resp. gene names
    names(l) <- c("foragers","nurses")
    
    # set the working directory to the pre-existing one
    setwd(current.dir)
    
    # return the list with all the plots and the data frame containing   
    return(l);
    
  }



# Stacked plots II --------------------------------------------------------

# Function name: stacked.zplot
# USage: to plot multiple diel gene expression (of Z scores) stacked on top of each other;
# Application: could be used to report similar diel expression patterns for multiple genes

stacked.logplot <- 
  function(gene_names, caste = "both", lwd=1.5, alpha=0.75, bg.alpha = 0.1) {
    
    # load the required libraries
    library(tidyverse)
    library(gridExtra)
    library(ggplot2)
    
    # load the core datasets
    load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_core_datasets.RData")
    
    #save the current directory 
    current.dir <- getwd()
    #change the directory to the folder where we have our zscores
    #setwd("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data")
    # save the list of gene names in a character vector
    g <- gene_names
    
    # if (log=F) {
    # # Read the zscores for caste
    # # Read the zscore file
    # # foragers_zscores <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data/TC5_zscores_foragers.csv", header = T, stringsAsFactors = F)
    # # nurses_zscores <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data/TC5_zscores_nurses.csv", header = T, stringsAsFactors = F)
    # foragers_zscores <- cflo.zscores.for
    # nurses_zscores <- cflo.zscores.nur
    # }

    # else {
      foragers_zscores <-  cflo.annots.exp %>% dplyr::select(gene_name, X2F:X24F)
      foragers_zscores[-1] <- log2(foragers_zscores[-1] + 1)
      nurses_zscores <-  cflo.annots.exp %>% dplyr::select(gene_name, X2N:X24N)
      nurses_zscores[-1] <- log2(nurses_zscores[-1] + 1)
    # }
    
    # Read the annotation file for description of the genes
    #all_genes <- read.csv("~/R-scripts_zombie_ant_lab/Enrichment_analysis/Ants/cflo_annotations_expression_v2.csv", header = T, stringsAsFactors = FALSE)
    all_genes <- cflo.annots.exp
    all_genes_annots <- all_genes %>% 
      dplyr::select(gene_name, annot = old_annotation)
    
    # Transforming the data to be able to use ggplot
    # Let's try the logic on a dummy subset
    dummy.for <- foragers_zscores %>%
      filter(gene_name %in% g) %>% 
      gather(ZT, zscore, -1) %>% 
      arrange(match(gene_name, g)) %>% 
      mutate(caste = "for") %>% 
      mutate(ZT = readr::parse_number(ZT)) %>% 
      left_join(all_genes_annots, by = "gene_name")
    
    dummy.nur <- nurses_zscores %>%
      filter(gene_name %in% g) %>% 
      gather(ZT, zscore, -1) %>% 
      arrange(match(gene_name, g)) %>% 
      mutate(caste = "nur") %>% 
      mutate(ZT = readr::parse_number(ZT)) %>% 
      left_join(all_genes_annots, by = "gene_name")
    
    # make an if else statement to modify the dataset as per the "caste" parameter
    if (caste == "for" | caste == "forager" | caste == "foragers" | caste == "f"){
      dummy <- dummy.for
      col.scheme <- c("#D91A2A")
    } 
    else if (caste == "nur" | caste == "nurse" | caste == "nurses" | caste == "n") {
      dummy <- dummy.nur
      col.scheme <- c("#2182BF")
    } else if (caste == "all" | caste == "both" | caste == "a" | caste == "b") {
      dummy <- rbind(dummy.for, dummy.nur)
      col.scheme <- c("#D91A2A","#2182BF")
    } else {
      print("Invalid value for caste. Use one of the following options: for, nur, all.")
      # stop the function and give the user an error.
      stop();
    }
    
    # make the gene_name and caste columns as factors
    dummy[[1]] <- as.factor(dummy[[1]]) # gene name
    dummy[[4]] <- as.factor(dummy[[4]]) # caste
    dummy[[5]] <- as.factor(dummy[[5]]) # annotation column
    
    dummy.summary <- 
      dummy %>% 
      group_by(caste, ZT) %>% 
      summarise(mean_zscore = mean(zscore))
    
    # Initialize a list to save the plots
    l <- list()
    # Let's plot
    library(ggplot2)
    pd <- position_dodge(0.2)
    l <- lapply(sort(unique(dummy[[4]])), function(i) {
      
      # Define the dataset
      ggplot() + 
        
        # # indicate light-dark phase
        # geom_rect(aes(xmin = 11.5, xmax = 23.5, ymin = -Inf, ymax = Inf),    
        #           fill = "lightgrey", alpha = 0.02, color=NA) +
        
        # Plot the individual gene expressions (zscores)
        geom_line(data=dummy[dummy$caste==i,], 
                  aes(x=as.numeric(as.character(ZT)), y=zscore, group = gene_name),
                  size=0.5, alpha=bg.alpha, 
                  col = ifelse(i=="for",
                               c(col.scheme[[1]], col.scheme[[2]]),
                               c(col.scheme[[2]], col.scheme[[1]]))) +
        
        # Plot the mean gene expression for the given gene list
        geom_line(data = dummy.summary[dummy.summary$caste == i,],
                  aes(x=as.numeric(as.character(ZT)), y=mean_zscore), 
                  col = "#3C3C40",
                  size = lwd, alpha = alpha) +
        
        
        # Set the theme
        xlab("") +
        ylab("") +
        theme_Publication() +
        scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
        scale_y_continuous(limits = c(min(dummy[[3]]),max(dummy[[3]]))) + 
        theme(text = element_text(size = 20, colour = 'black'),
              legend.position = "none",
              axis.title.x=element_blank(),
              axis.text.x=element_blank()) +
        # set transparency
        theme( 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA)) 
      # show only a subset of the y-axis tick labels
      # theme(axis.text.y=element_text(color=c("transparent","black","transparent",
      #                                        "transparent","black","transparent",
      #                                        "transparent","transparent","transparent",
      #                                        "black")))
    })
    
    # Let's name the plots with their resp. gene names
    names(l) <- c("foragers","nurses")
    
    # set the working directory to the pre-existing one
    setwd(current.dir)
    
    # return the list with all the plots and the data frame containing   
    return(l);
    
  }

