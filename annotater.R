#!/usr/bin/env Rscript

options(warn = -1)

# Load libraries

suppressPackageStartupMessages({
library (ggplot2)
library (dplyr)
library (optparse)
library (ggrepel)
library (tidyr)
})



#-------------------------------------------------------
# Function: doCoverage. Creates coverage/dentity plots
#-------------------------------------------------------


doCoverage <- function (x)
{

  dfAgo.i <- dfAgo %>% filter (name == x)
 
  
  nflank <- (nchar  (dfAgo.i$seq) - (abs(dfAgo.i$end - dfAgo.i$start))) / 2  
    
  all_pos <- data.frame(pos = (dfAgo.i$start-nflank):(dfAgo.i$end+nflank-1))
  
  # make sure all positions are represented in the plot
  
  dfCov.i <- merge (dfCov %>% dplyr::filter (id==x) %>% arrange (length), all_pos, by="pos", all.y=TRUE)
        
  # convert NA to 0
  
  dfCov.i$reads[is.na(dfCov.i$reads)] <- 0
  
  # arrange data.frame by read-length
  
  dfCov.i <-
    dfCov.i %>%
    group_by (pos) %>%
    arrange (length)
  
  
  label <- paste (dfAgo.i$host_name, " (", dfAgo.i$X.chrom, ":", dfAgo.i$start, "-", dfAgo.i$end, ")", sep="")
  
  # ggplot
  
 gg <- 
    ggplot () + 
    geom_bar (data=dfCov.i, aes (x=pos, y=reads, fill=length), stat="identity") + 
    scale_fill_gradient2(low = color_low, mid=color_mid, high=color_high, midpoint=30) +     
    labs (x="", y="Total Reads", fill="Read length") +         
    myTheme +
    theme (axis.text.x = element_blank(), axis.ticks.x = element_blank ())
  
  
  if (dfAgo.i$strand == "-") {
  
    gg <- gg + scale_x_reverse ()
  
  }
  
  # Add sequence below plot
  
  yrange <- ggplot_build(gg)$layout$panel_ranges[[1]]$y.range
  y <- -(yrange[2]) * 0.015

  if (dfAgo.i$strand == "+") {  
    gg <- gg + geom_text (data=NULL, aes (x=-Inf, y=Inf, label=label), hjust=-0.02, vjust=1.2)
    dfSeq <- data.frame (seq = unlist(strsplit (dfAgo.i$seq, "")), pos = all_pos, y=y)
  } else {
    gg <- gg + geom_text (data=NULL, aes (x=Inf, y=Inf, label=label), hjust=-0.02, vjust=1.2)
    dfSeq <- data.frame (seq = rev(unlist(strsplit (dfAgo.i$seq, ""))), pos = all_pos, y=y)
  }
  
  gg <- gg + geom_text (data=dfSeq, aes (x=pos, y=y, label=seq), size=2)

  print (gg)
  
}


#-------------------------------------------------------
# Function: doExpression. Creates expression plots
#-------------------------------------------------------


doExpression <- function (x)
{

   columns = "RPMM"

   cn <- colnames (dfAgo)
  
   cn.exp <- cn[grep (paste (columns, "_", sep=""), cn)]
  
   dfAgo.i <- dfAgo %>% filter (name == x)
   
   dfAgo.m <- 
     dfAgo.i %>%
     select (c(cn[1:7], cn.exp)) %>%
     gather (Sample, RPMM, cn.exp)
   
   dfAgo.m$label <- gsub(paste (columns, '_(.*?)(_).*', sep=""), '\\1', gsub ("\\.", "_", dfAgo.m$Sample))
  
   locus <- paste (dfAgo.i$host_name, " (", dfAgo.i$X.chrom, ":", dfAgo.i$start, "-", dfAgo.i$end, ")", sep="")
   
   gg <-
     ggplot (dfAgo.m, aes (x=label, y=RPMM)) +
     geom_bar (stat="identity") +     
     geom_text (data=NULL, aes (x=-Inf, y=Inf, label=locus), hjust=-0.02, vjust=1.2) + 
     labs (x="Datasets", y="RPMM") +      
     myTheme
    
  
   print (gg)
   
}



myTheme <- theme_bw() + 
           theme (axis.title = element_text(size = 18, color = "black"),
                  axis.text.y = element_text(size = 16, color = "black"),
                  axis.text.x = element_text(size = 14, color = "black"),
                  legend.title = element_text(size = 17, face = NULL, color = "black"),
                  legend.text = element_text(size = 16, face = NULL, color = "black"),
                  legend.background = element_rect(fill = "transparent", colour = "transparent"),
                  panel.background = element_rect(fill = "white", color="grey", size=1))



#-------------------------------------------------------
# Options
#-------------------------------------------------------

version <- "1.0.0" 
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="Input agotron data file", metavar="file"),
  make_option(c("-c", "--coverage"), type="character", default="coverage.txt", help="Input coverage file", metavar="file"),
  make_option(c("-p", "--prefix"), type="character", default="agotron", help="output prefix [default= %default]", metavar="string"),  
  make_option(c("-m", "--min_median"), type="integer", default=30, help="Lowest median read length to be classified as agotrom [default= %default]", metavar="int"),
  make_option(c("-e", "--min_homogeneity"), type="numeric", default=0.7, help="Fraction of homogenous 5' ends [default= %default]", metavar="float"),
  make_option(c("-d", "--min_distance"), type="integer", default=1, help="Minimun distance from prevalent 5'end to 5'ss [default= %default]", metavar="int")
); 

opt_parser = OptionParser(option_list=option_list,                           
                          description = "Annotate and visualize agotrons from analyzer.py output");
opt = parse_args(opt_parser);

color_high <- "#4f94d8"
color_low <- "#dc4855"
color_mid <- "grey"


#-------------------------------------------------------
# Load data
#-------------------------------------------------------

dfAgo <- read.table (file("stdin"), sep="\t", header=T, stringsAsFactors=F, comment.char="!")

if (!is.null (opt$coverage) && file.exists(opt$coverage)) {
   dfCov <- read.table (opt$coverage, sep="\t", header=T, stringsAsFactors=F, comment.char="!")
}



dfAgo$p5distance <- ifelse (dfAgo$strand == "+", abs(dfAgo$p5position - dfAgo$start),  abs(dfAgo$p5position - dfAgo$end))
dfAgo$p5agotron <- dfAgo$p5homogeneity * ifelse (dfAgo$p5distance <= opt$min_distance, 1, 0) #

dfAgo$isAgotron <- (dfAgo$median > opt$min_median & 
                    dfAgo$p5homogeneity > opt$min_homogeneity & 
                    dfAgo$p5distance <= opt$min_distance)
                    

                    

#-------------------------------------------------------
# Print agotrons
#-------------------------------------------------------
dfAgotrons <- subset (dfAgo, dfAgo$isAgotron == T)

print (paste ("Found", nrow (dfAgotrons), "agotrons"))

filename <- paste (opt$prefix, "filtered.bed", sep="_")

print (paste ("Writing...", filename))

cn <- colnames (dfAgotrons)  
output_columns <- c(cn[c(1:3,7,5,6)], cn[grep (paste ("RPMM", "_", sep=""), cn)], cn[grep (paste ("reads", "_", sep=""), cn)]) 
         
         
write.table (dfAgotrons[order (dfAgotrons$score, decreasing=T),output_columns], filename, sep="\t", row.names = F, quote=F)

                                      
                    

#-------------------------------------------------------
# Scatter plot
#-------------------------------------------------------

filename <- paste (opt$prefix, "scatter.pdf", sep="_")

print (paste ("Plotting...", filename))

nlabels <- 10

# label top20 expressed introns

toLabel <-  dfAgo %>% arrange (-score) %>%  "["(.,1:nlabels,"name") 


# OR label top20 expressed agotrons

toLabel <-  dfAgo %>% arrange (-isAgotron, -score) %>%  "["(.,1:nlabels,"name") 




pdf (file=filename, width=10)

ggplot (dfAgo, aes (x=score, y=median)) + 
  geom_errorbar(aes(ymin=percentile_5th, ymax=percentile_95th), width=.05, col="grey") +
  geom_point (size=3, aes(color = p5agotron)) + 
  scale_color_gradient2(low = color_low, mid=color_mid, high=color_high, midpoint=0.5, breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%")) + 
  geom_text_repel(aes(label = ifelse(name %in% toLabel, as.character(host_name),'')), box.padding = unit(0.45, "lines")) +
  scale_x_log10 () + 
  labs(x="Total reads", y="Median read length", colour = "5'end@SD") +  
  myTheme
  

  
invisible (dev.off())




#-------------------------------------------------------
# Additional plots
#-------------------------------------------------------

nplot <- 10

# Plot top10 expressed agotrons
toPlot <-  dfAgo %>% arrange (-isAgotron, -score) %>%  "["(.,1:nplot,"name") 


# Or top10 expressed loci
toPlot <-  dfAgo %>% arrange (-score) %>%  "["(.,1:nplot,"name") 



#-------------------------------------------------------
# Expression plot
#-------------------------------------------------------


filename <- paste (opt$prefix, "expression.pdf", sep="_")

print (paste ("Plotting...", filename))

pdf (file=filename, width=10, height=6)

#mapply (doExpression, toPlot, SIMPLIFY = TRUE)

invisible (sapply (toPlot, doExpression))
invisible (dev.off())

#-------------------------------------------------------
# Coverage/Density plot
#-------------------------------------------------------

if (!is.null (opt$coverage)) {

    
    filename <- paste (opt$prefix, "coverage.pdf", sep="_")

    print (paste ("Plotting...", filename))
            
    pdf (file=filename, width=10, height=6)
    #a <- mapply (doCoverage, toPlot)
    invisible(sapply (toPlot, doCoverage))
    invisible (dev.off() )

}

